#!/usr/bin/env python3
"""
SOLID Principles Analyzer
Analyzes Python code for SOLID principle violations.
"""

import ast
import sys
from dataclasses import dataclass
from typing import List, Dict, Any
from pathlib import Path


@dataclass
class Violation:
    """Represents a SOLID principle violation."""
    principle: str  # S, O, L, I, or D
    line: int
    class_name: str
    message: str
    severity: str  # 'high', 'medium', 'low'
    suggestion: str


class SOLIDAnalyzer(ast.NodeVisitor):
    """AST visitor that detects SOLID violations."""

    def __init__(self, source_code: str):
        self.source_code = source_code
        self.violations: List[Violation] = []
        self.current_class = None
        self.lines = source_code.split('\n')

    def analyze(self) -> List[Violation]:
        """Run the analysis."""
        tree = ast.parse(self.source_code)
        self.visit(tree)
        return self.violations

    def visit_ClassDef(self, node: ast.ClassDef):
        """Analyze class definitions."""
        self.current_class = node.name

        # Check Single Responsibility Principle
        self._check_single_responsibility(node)

        # Check Interface Segregation Principle
        self._check_interface_segregation(node)

        # Check Liskov Substitution Principle
        self._check_liskov_substitution(node)

        # Check Dependency Inversion Principle
        self._check_dependency_inversion(node)

        self.generic_visit(node)
        self.current_class = None

    def visit_FunctionDef(self, node: ast.FunctionDef):
        """Analyze function definitions."""
        if self.current_class:
            # Check for Open/Closed violations in methods
            self._check_open_closed(node)

        self.generic_visit(node)

    def _check_single_responsibility(self, node: ast.ClassDef):
        """Check for Single Responsibility violations."""
        # Count lines in class
        class_lines = node.end_lineno - node.lineno

        if class_lines > 100:
            self.violations.append(Violation(
                principle='S',
                line=node.lineno,
                class_name=node.name,
                message=f'Class has {class_lines} lines (>100). Likely has multiple responsibilities.',
                severity='high',
                suggestion='Apply Extract Class pattern to separate concerns.'
            ))

        # Count number of methods (excluding __init__, __str__, etc.)
        methods = [n for n in node.body if isinstance(n, ast.FunctionDef)]
        public_methods = [m for m in methods if not m.name.startswith('_')]

        if len(public_methods) > 10:
            self.violations.append(Violation(
                principle='S',
                line=node.lineno,
                class_name=node.name,
                message=f'Class has {len(public_methods)} public methods. Too many responsibilities.',
                severity='medium',
                suggestion='Split into multiple focused classes, each with a single responsibility.'
            ))

        # Check for mixed concerns (DB + UI + Business Logic)
        method_names = [m.name for m in methods]
        has_db = any('save' in m or 'load' in m or 'query' in m for m in method_names)
        has_ui = any('render' in m or 'display' in m or 'show' in m for m in method_names)
        has_validation = any('validate' in m or 'check' in m for m in method_names)

        mixed_concerns = sum([has_db, has_ui, has_validation])
        if mixed_concerns >= 2:
            self.violations.append(Violation(
                principle='S',
                line=node.lineno,
                class_name=node.name,
                message='Class mixes multiple concerns (data access, UI, validation).',
                severity='high',
                suggestion='Separate into Repository (data), View (UI), and Validator (logic) classes.'
            ))

    def _check_open_closed(self, node: ast.FunctionDef):
        """Check for Open/Closed violations."""
        # Look for long if/elif chains that switch on type
        for child in ast.walk(node):
            if isinstance(child, ast.If):
                elif_count = 0
                current = child

                while hasattr(current, 'orelse') and current.orelse:
                    if isinstance(current.orelse[0], ast.If):
                        elif_count += 1
                        current = current.orelse[0]
                    else:
                        break

                if elif_count >= 3:
                    self.violations.append(Violation(
                        principle='O',
                        line=child.lineno,
                        class_name=self.current_class or 'module',
                        message=f'Long if/elif chain ({elif_count + 1} branches). Closed for extension.',
                        severity='high',
                        suggestion='Replace with Strategy pattern or polymorphism.'
                    ))

    def _check_liskov_substitution(self, node: ast.ClassDef):
        """Check for Liskov Substitution violations."""
        methods = [n for n in node.body if isinstance(n, ast.FunctionDef)]

        for method in methods:
            # Check for empty overrides or NotImplementedError
            if len(method.body) == 1:
                stmt = method.body[0]

                # Check for "pass"
                if isinstance(stmt, ast.Pass):
                    self.violations.append(Violation(
                        principle='L',
                        line=method.lineno,
                        class_name=node.name,
                        message=f'Method {method.name} is empty override. Breaks LSP.',
                        severity='high',
                        suggestion='Redesign inheritance hierarchy. Child should fully implement parent behavior.'
                    ))

                # Check for "raise NotImplementedError"
                if isinstance(stmt, ast.Raise):
                    if isinstance(stmt.exc, ast.Call):
                        if isinstance(stmt.exc.func, ast.Name):
                            if stmt.exc.func.id == 'NotImplementedError':
                                self.violations.append(Violation(
                                    principle='L',
                                    line=method.lineno,
                                    class_name=node.name,
                                    message=f'Method {method.name} raises NotImplementedError. Violates LSP.',
                                    severity='high',
                                    suggestion='Use abstract base class or redesign to avoid forcing empty implementations.'
                                ))

    def _check_interface_segregation(self, node: ast.ClassDef):
        """Check for Interface Segregation violations."""
        # If class inherits from multiple bases and has many methods
        if len(node.bases) > 0:
            methods = [n for n in node.body if isinstance(n, ast.FunctionDef)]

            if len(methods) > 15:
                self.violations.append(Violation(
                    principle='I',
                    line=node.lineno,
                    class_name=node.name,
                    message=f'Class/Interface has {len(methods)} methods. Too fat.',
                    severity='medium',
                    suggestion='Split into multiple focused interfaces/protocols.'
                ))

    def _check_dependency_inversion(self, node: ast.ClassDef):
        """Check for Dependency Inversion violations."""
        # Look for direct instantiation in __init__
        for method in node.body:
            if isinstance(method, ast.FunctionDef) and method.name == '__init__':
                for stmt in ast.walk(method):
                    if isinstance(stmt, ast.Assign):
                        for target in stmt.targets:
                            if isinstance(target, ast.Attribute):
                                # Check if assigning to self.something = SomeClass()
                                if isinstance(stmt.value, ast.Call):
                                    if isinstance(stmt.value.func, ast.Name):
                                        # This is direct instantiation
                                        class_being_created = stmt.value.func.id

                                        # Common violation patterns
                                        if any(keyword in class_being_created.lower()
                                               for keyword in ['database', 'repository', 'service', 'client', 'api']):
                                            self.violations.append(Violation(
                                                principle='D',
                                                line=stmt.lineno,
                                                class_name=node.name,
                                                message=f'Direct instantiation of {class_being_created} in __init__. Depends on concrete class.',
                                                severity='high',
                                                suggestion='Inject dependency via __init__ parameter with Protocol type.'
                                            ))


def analyze_file(file_path: str) -> List[Violation]:
    """Analyze a Python file for SOLID violations."""
    with open(file_path, 'r', encoding='utf-8') as f:
        source_code = f.read()

    analyzer = SOLIDAnalyzer(source_code)
    return analyzer.analyze()


def print_violations(violations: List[Violation], file_path: str = None):
    """Print violations in a readable format."""
    if not violations:
        print("✅ No SOLID violations detected!")
        return

    print(f"\n{'='*80}")
    if file_path:
        print(f"SOLID Analysis: {file_path}")
    print(f"{'='*80}\n")

    # Group by principle
    by_principle = {}
    for v in violations:
        if v.principle not in by_principle:
            by_principle[v.principle] = []
        by_principle[v.principle].append(v)

    principle_names = {
        'S': 'Single Responsibility',
        'O': 'Open/Closed',
        'L': 'Liskov Substitution',
        'I': 'Interface Segregation',
        'D': 'Dependency Inversion'
    }

    for principle in ['S', 'O', 'L', 'I', 'D']:
        if principle in by_principle:
            print(f"\n[{principle}] {principle_names[principle]} Violations:")
            print("-" * 80)

            for v in by_principle[principle]:
                severity_icon = {
                    'high': '🔴',
                    'medium': '🟡',
                    'low': '🟢'
                }[v.severity]

                print(f"\n{severity_icon} Line {v.line} in {v.class_name}")
                print(f"   Problem: {v.message}")
                print(f"   Fix: {v.suggestion}")

    print(f"\n{'='*80}")
    print(f"Total violations: {len(violations)}")
    print(f"{'='*80}\n")


def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print("Usage: python solid_analyzer.py <file.py>")
        sys.exit(1)

    file_path = sys.argv[1]

    if not Path(file_path).exists():
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)

    violations = analyze_file(file_path)
    print_violations(violations, file_path)

    # Exit with error code if violations found
    sys.exit(1 if violations else 0)


if __name__ == '__main__':
    main()
