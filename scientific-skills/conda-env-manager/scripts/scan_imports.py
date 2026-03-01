#!/usr/bin/env python3
"""
Project Import Scanner

Scans a Python project directory for all third-party imports using AST parsing.
Outputs a JSON report with:
  - required: hard imports
  - optional: imports inside try/except blocks
  - missing: imports not installed in the target conda environment
"""

import ast
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, Set, Tuple

# Python standard library modules (3.10-3.13)
STDLIB_MODULES = {
    "__future__", "_thread", "abc", "aifc", "argparse", "array", "ast",
    "asynchat", "asyncio", "asyncore", "atexit", "audioop", "base64",
    "bdb", "binascii", "binhex", "bisect", "builtins", "bz2", "calendar",
    "cgi", "cgitb", "chunk", "cmath", "cmd", "code", "codecs", "codeop",
    "collections", "colorsys", "compileall", "concurrent", "configparser",
    "contextlib", "contextvars", "copy", "copyreg", "cProfile", "crypt",
    "csv", "ctypes", "curses", "dataclasses", "datetime", "dbm", "decimal",
    "difflib", "dis", "distutils", "doctest", "email", "encodings",
    "enum", "errno", "faulthandler", "fcntl", "filecmp", "fileinput",
    "fnmatch", "fractions", "ftplib", "functools", "gc", "getopt",
    "getpass", "gettext", "glob", "graphlib", "grp", "gzip", "hashlib",
    "heapq", "hmac", "html", "http", "idlelib", "imaplib", "imghdr",
    "imp", "importlib", "inspect", "io", "ipaddress", "itertools", "json",
    "keyword", "lib2to3", "linecache", "locale", "logging", "lzma",
    "mailbox", "mailcap", "marshal", "math", "mimetypes", "mmap",
    "modulefinder", "multiprocessing", "netrc", "nis", "nntplib",
    "numbers", "operator", "optparse", "os", "ossaudiodev", "pathlib",
    "pdb", "pickle", "pickletools", "pipes", "pkgutil", "platform",
    "plistlib", "poplib", "posix", "posixpath", "pprint", "profile",
    "pstats", "pty", "pwd", "py_compile", "pyclbr", "pydoc", "queue",
    "quopri", "random", "re", "readline", "reprlib", "resource", "rlcompleter",
    "runpy", "sched", "secrets", "select", "selectors", "shelve", "shlex",
    "shutil", "signal", "site", "smtpd", "smtplib", "sndhdr", "socket",
    "socketserver", "spwd", "sqlite3", "sre_compile", "sre_constants",
    "sre_parse", "ssl", "stat", "statistics", "string", "stringprep",
    "struct", "subprocess", "sunau", "symtable", "sys", "sysconfig",
    "syslog", "tabnanny", "tarfile", "telnetlib", "tempfile", "termios",
    "test", "textwrap", "threading", "time", "timeit", "tkinter", "token",
    "tokenize", "tomllib", "trace", "traceback", "tracemalloc", "tty",
    "turtle", "turtledemo", "types", "typing", "unicodedata", "unittest",
    "urllib", "uu", "uuid", "venv", "warnings", "wave", "weakref",
    "webbrowser", "winreg", "winsound", "wsgiref", "xdrlib", "xml",
    "xmlrpc", "zipapp", "zipfile", "zipimport", "zlib", "zoneinfo",
    # Common sub-modules that appear as top-level
    "collections.abc", "concurrent.futures", "email.mime", "http.client",
    "http.server", "importlib.metadata", "multiprocessing.pool",
    "os.path", "typing_extensions", "unittest.mock", "urllib.parse",
    "urllib.request", "xml.etree", "xml.etree.ElementTree",
}

# Import name → pip/conda package name mapping
IMPORT_TO_PACKAGE = {
    "PIL": "Pillow",
    "cv2": "opencv-python",
    "sklearn": "scikit-learn",
    "skopt": "scikit-optimize",
    "Bio": "biopython",
    "rdkit": "rdkit",  # conda-forge only
    "yaml": "pyyaml",
    "bs4": "beautifulsoup4",
    "attr": "attrs",
    "wx": "wxPython",
    "gi": "PyGObject",
    "Crypto": "pycryptodome",
    "serial": "pyserial",
    "usb": "pyusb",
    "zmq": "pyzmq",
    "dotenv": "python-dotenv",
    "jose": "python-jose",
    "magic": "python-magic",
    "dateutil": "python-dateutil",
    "lxml": "lxml",
    "mpl_toolkits": "matplotlib",
    "google": "google-api-python-client",
    "umap": "umap-learn",
    "adjustText": "adjustText",
    "CoolProp": "CoolProp",
}

# Packages that should be installed via conda-forge (not pip)
CONDA_FORGE_ONLY = {
    "rdkit", "openmm", "openff", "openmmforcefields", "pdbfixer",
    "pymol", "numpoly",
}


class ImportVisitor(ast.NodeVisitor):
    """AST visitor that extracts import statements."""

    def __init__(self):
        self.required: Set[str] = set()
        self.optional: Set[str] = set()
        self._in_try = False

    def visit_Try(self, node: ast.Try):
        # Check if any except handler catches ImportError
        catches_import = any(
            h.type and (
                (isinstance(h.type, ast.Name) and h.type.id in ("ImportError", "ModuleNotFoundError"))
                or (isinstance(h.type, ast.Tuple) and any(
                    isinstance(e, ast.Name) and e.id in ("ImportError", "ModuleNotFoundError")
                    for e in h.type.elts
                ))
            )
            for h in node.handlers
        )

        if catches_import:
            old = self._in_try
            self._in_try = True
            for child in node.body:
                self.visit(child)
            self._in_try = old
            for handler in node.handlers:
                for child in handler.body:
                    self.visit(child)
        else:
            self.generic_visit(node)

    def visit_Import(self, node: ast.Import):
        for alias in node.names:
            top = alias.name.split(".")[0]
            if self._in_try:
                self.optional.add(top)
            else:
                self.required.add(top)

    def visit_ImportFrom(self, node: ast.ImportFrom):
        if node.module:
            top = node.module.split(".")[0]
            if self._in_try:
                self.optional.add(top)
            else:
                self.required.add(top)


def scan_project(project_dir: str) -> Dict:
    """Scan all .py files in a project for imports."""
    root = Path(project_dir)
    all_required: Set[str] = set()
    all_optional: Set[str] = set()
    file_count = 0
    errors = []

    # Directories to skip
    skip_dirs = {
        ".git", "__pycache__", ".tox", ".eggs", "node_modules",
        ".venv", "venv", ".idea", ".vscode", ".claude",
    }

    for py_file in root.rglob("*.py"):
        # Skip excluded directories
        if any(part in skip_dirs for part in py_file.parts):
            continue
        try:
            source = py_file.read_text(encoding="utf-8", errors="ignore")
            tree = ast.parse(source, filename=str(py_file))
            visitor = ImportVisitor()
            visitor.visit(tree)
            all_required.update(visitor.required)
            all_optional.update(visitor.optional)
            file_count += 1
        except SyntaxError:
            errors.append(str(py_file.relative_to(root)))

    # Remove stdlib
    third_party_required = {m for m in all_required if m not in STDLIB_MODULES}
    third_party_optional = {m for m in all_optional if m not in STDLIB_MODULES}

    # Remove project-internal modules
    # 1) Top-level packages and .py files
    project_modules = set()
    for item in root.iterdir():
        if item.is_dir() and (item / "__init__.py").exists():
            project_modules.add(item.name)
        elif item.suffix == ".py" and item.stem != "__init__":
            project_modules.add(item.stem)

    # 2) All .py file stems and subdirectory names in the project
    #    (catches internal modules like _boiler_turbogenerator, evaluation, units)
    for py_file in root.rglob("*.py"):
        if any(part in skip_dirs for part in py_file.parts):
            continue
        stem = py_file.stem
        if stem != "__init__":
            project_modules.add(stem)

    # 3) All subdirectory names (even without __init__.py, for namespace packages)
    for item in root.rglob("*"):
        if item.is_dir() and not any(part in skip_dirs for part in item.parts):
            project_modules.add(item.name)

    # 4) Private modules starting with _ are almost always internal
    private_modules = {m for m in (third_party_required | third_party_optional)
                       if m.startswith("_") and m not in STDLIB_MODULES}
    project_modules.update(private_modules)

    third_party_required -= project_modules
    third_party_optional -= project_modules

    # Optional imports that are also required elsewhere → required
    third_party_optional -= third_party_required

    # Map to package names
    def map_packages(imports: Set[str]) -> Dict[str, str]:
        result = {}
        for imp in sorted(imports):
            pkg = IMPORT_TO_PACKAGE.get(imp, imp)
            result[imp] = pkg
        return result

    return {
        "project": str(root),
        "files_scanned": file_count,
        "parse_errors": errors,
        "required": map_packages(third_party_required),
        "optional": map_packages(third_party_optional),
        "conda_forge_only": sorted(
            imp for imp in (third_party_required | third_party_optional)
            if IMPORT_TO_PACKAGE.get(imp, imp) in CONDA_FORGE_ONLY
            or imp in CONDA_FORGE_ONLY
        ),
    }


def check_installed(env_name: str, imports: Dict[str, str]) -> Tuple[Dict, Dict]:
    """Check which packages are installed in a conda environment."""
    try:
        result = subprocess.run(
            ["conda", "list", "-n", env_name, "--json"],
            capture_output=True, text=True, timeout=30
        )
        installed_raw = json.loads(result.stdout) if result.returncode == 0 else []
    except Exception:
        installed_raw = []

    installed_names = {pkg["name"].lower() for pkg in installed_raw}

    found = {}
    missing = {}
    for imp, pkg in imports.items():
        # Check both import name and package name (lowercased)
        if imp.lower() in installed_names or pkg.lower() in installed_names:
            found[imp] = pkg
        else:
            missing[imp] = pkg

    return found, missing


def main():
    if len(sys.argv) < 2:
        print("Usage: scan_imports.py <project_dir> [conda_env_name]")
        sys.exit(1)

    project_dir = sys.argv[1]
    env_name = sys.argv[2] if len(sys.argv) > 2 else None

    report = scan_project(project_dir)

    if env_name:
        all_imports = {**report["required"], **report["optional"]}
        found, missing = check_installed(env_name, all_imports)
        report["env_name"] = env_name
        report["installed"] = found
        report["missing"] = missing

    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
