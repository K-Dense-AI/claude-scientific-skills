#!/usr/bin/env python3
"""
Project Analyzer - 프로젝트 분석 템플릿

이 스크립트는 프로젝트 구조와 코드 품질을 분석하는 템플릿입니다.
실제 구현은 Claude Code의 도구(Glob, Grep, Read)를 사용합니다.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Dict
import re


@dataclass
class FileAnalysis:
    """파일 분석 결과"""
    path: str
    lines: int
    classes: List['ClassAnalysis']
    functions: int
    imports: List[str]


@dataclass
class ClassAnalysis:
    """클래스 분석 결과"""
    name: str
    lines: int
    methods: int
    file_path: str


@dataclass
class ProjectReport:
    """프로젝트 분석 리포트"""
    total_files: int
    total_lines: int
    critical_files: List[FileAnalysis]  # > 800줄
    moderate_files: List[FileAnalysis]  # > 500줄
    critical_classes: List[ClassAnalysis]  # > 450줄
    moderate_classes: List[ClassAnalysis]  # > 300줄
    duplicates: List[Tuple[str, str]]  # 중복 파일명
    missing_init: List[str]  # __init__.py 없는 디렉토리
    solid_violations: Dict[str, List[str]]


class ProjectAnalyzer:
    """프로젝트 분석기"""

    def __init__(self, project_root: str):
        self.root = Path(project_root)

    def analyze(self) -> ProjectReport:
        """전체 프로젝트 분석"""
        files = self._find_python_files()
        file_analyses = [self._analyze_file(f) for f in files]

        return ProjectReport(
            total_files=len(files),
            total_lines=sum(f.lines for f in file_analyses),
            critical_files=self._filter_critical_files(file_analyses),
            moderate_files=self._filter_moderate_files(file_analyses),
            critical_classes=self._filter_critical_classes(file_analyses),
            moderate_classes=self._filter_moderate_classes(file_analyses),
            duplicates=self._find_duplicates(files),
            missing_init=self._find_missing_init(),
            solid_violations=self._check_solid_violations(file_analyses),
        )

    def _find_python_files(self) -> List[Path]:
        """Python 파일 찾기"""
        # 실제 구현: Glob tool 사용
        # Glob(pattern="**/*.py", path=self.root)
        return list(self.root.rglob("*.py"))

    def _analyze_file(self, file_path: Path) -> FileAnalysis:
        """개별 파일 분석"""
        # 실제 구현: Read tool 사용
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
            lines = content.count('\n') + 1

        classes = self._extract_classes(file_path, content)
        functions = len(re.findall(r'^def \w+\(', content, re.MULTILINE))
        imports = self._extract_imports(content)

        return FileAnalysis(
            path=str(file_path),
            lines=lines,
            classes=classes,
            functions=functions,
            imports=imports,
        )

    def _extract_classes(self, file_path: Path, content: str) -> List[ClassAnalysis]:
        """클래스 추출 및 분석"""
        classes = []
        class_pattern = r'^class (\w+).*?:'
        matches = re.finditer(class_pattern, content, re.MULTILINE)

        for match in matches:
            class_name = match.group(1)
            # 클래스의 줄 수 계산 (단순화된 버전)
            class_lines = self._count_class_lines(content, match.start())
            method_count = self._count_methods(content, match.start(), class_lines)

            classes.append(ClassAnalysis(
                name=class_name,
                lines=class_lines,
                methods=method_count,
                file_path=str(file_path),
            ))

        return classes

    def _count_class_lines(self, content: str, start_pos: int) -> int:
        """클래스의 줄 수 계산 (단순화)"""
        # 실제로는 AST 파싱이 더 정확함
        lines = content[start_pos:].split('\n')
        indent_level = len(lines[0]) - len(lines[0].lstrip())

        count = 0
        for line in lines[1:]:
            if line.strip() and not line.startswith(' ' * (indent_level + 1)):
                break
            count += 1

        return count

    def _count_methods(self, content: str, start_pos: int, class_lines: int) -> int:
        """클래스의 메서드 수 계산"""
        class_content = '\n'.join(content[start_pos:].split('\n')[:class_lines])
        return len(re.findall(r'^\s+def \w+\(', class_content, re.MULTILINE))

    def _extract_imports(self, content: str) -> List[str]:
        """import 문 추출"""
        imports = []
        import_pattern = r'^(?:from|import)\s+(\w+)'
        matches = re.finditer(import_pattern, content, re.MULTILINE)
        imports = [match.group(1) for match in matches]
        return imports

    def _filter_critical_files(self, files: List[FileAnalysis]) -> List[FileAnalysis]:
        """Critical 파일 필터링 (> 800줄)"""
        return [f for f in files if f.lines > 800]

    def _filter_moderate_files(self, files: List[FileAnalysis]) -> List[FileAnalysis]:
        """Moderate 파일 필터링 (> 500줄)"""
        return [f for f in files if 500 < f.lines <= 800]

    def _filter_critical_classes(self, files: List[FileAnalysis]) -> List[ClassAnalysis]:
        """Critical 클래스 필터링 (> 450줄)"""
        all_classes = [cls for f in files for cls in f.classes]
        return [cls for cls in all_classes if cls.lines > 450]

    def _filter_moderate_classes(self, files: List[FileAnalysis]) -> List[ClassAnalysis]:
        """Moderate 클래스 필터링 (> 300줄)"""
        all_classes = [cls for f in files for cls in f.classes]
        return [cls for cls in all_classes if 300 < cls.lines <= 450]

    def _find_duplicates(self, files: List[Path]) -> List[Tuple[str, str]]:
        """중복 파일명 찾기"""
        file_names = {}
        duplicates = []

        for file_path in files:
            name = file_path.name
            if name in file_names:
                duplicates.append((str(file_names[name]), str(file_path)))
            else:
                file_names[name] = file_path

        return duplicates

    def _find_missing_init(self) -> List[str]:
        """__init__.py 없는 디렉토리 찾기"""
        missing = []
        for directory in self.root.rglob('*'):
            if directory.is_dir():
                init_file = directory / '__init__.py'
                if not init_file.exists():
                    # Python 패키지 디렉토리인지 확인
                    if any(directory.glob('*.py')):
                        missing.append(str(directory))
        return missing

    def _check_solid_violations(self, files: List[FileAnalysis]) -> Dict[str, List[str]]:
        """SOLID 위반 사전 검사 (개요만)"""
        violations = {
            'god_classes': [],
            'manager_classes': [],
            'direct_dependencies': [],
            'isinstance_checks': [],
        }

        # 실제 구현: Grep tool 사용
        # 여기서는 간단한 패턴만 체크
        for file_analysis in files:
            # God 클래스: 30개 이상의 메서드
            for cls in file_analysis.classes:
                if cls.methods > 30:
                    violations['god_classes'].append(f"{cls.file_path}:{cls.name}")

                # Manager 패턴
                if 'Manager' in cls.name or 'Handler' in cls.name:
                    violations['manager_classes'].append(f"{cls.file_path}:{cls.name}")

            # 구체 클래스 직접 import (간단한 휴리스틱)
            for imp in file_analysis.imports:
                if 'Database' in imp or 'Connection' in imp:
                    violations['direct_dependencies'].append(f"{file_analysis.path}:{imp}")

        return violations

    def generate_report(self, report: ProjectReport) -> str:
        """마크다운 리포트 생성"""
        md = []
        md.append("# 프로젝트 분석 리포트\n")

        # 요약
        md.append("## 요약\n")
        md.append(f"- 총 파일: {report.total_files}개")
        md.append(f"- 총 줄 수: {report.total_lines:,}줄")
        md.append(f"- Critical 파일: {len(report.critical_files)}개")
        md.append(f"- Critical 클래스: {len(report.critical_classes)}개\n")

        # Critical 이슈
        md.append("## 🔴 Critical Issues\n")

        if report.critical_files:
            md.append("### 파일 크기 초과 (> 800줄)\n")
            for f in report.critical_files:
                md.append(f"- {f.path}: {f.lines}줄")
            md.append("")

        if report.critical_classes:
            md.append("### 클래스 크기 초과 (> 450줄)\n")
            for cls in report.critical_classes:
                md.append(f"- {cls.file_path}:{cls.name}: {cls.lines}줄 ({cls.methods}개 메서드)")
            md.append("")

        if report.missing_init:
            md.append("### __init__.py 누락\n")
            for d in report.missing_init:
                md.append(f"- {d}")
            md.append("")

        # Moderate 이슈
        md.append("## 🟡 Moderate Issues\n")

        if report.moderate_files:
            md.append("### 파일 크기 주의 (> 500줄)\n")
            for f in report.moderate_files:
                md.append(f"- {f.path}: {f.lines}줄")
            md.append("")

        if report.moderate_classes:
            md.append("### 클래스 크기 주의 (> 300줄)\n")
            for cls in report.moderate_classes:
                md.append(f"- {cls.file_path}:{cls.name}: {cls.lines}줄")
            md.append("")

        # SOLID 위반
        md.append("## ⚠️ SOLID 위반 가능성 (개요)\n")
        md.append("*상세 검증은 Phase 3에서 solid-principles 스킬이 수행합니다.*\n")

        for violation_type, items in report.solid_violations.items():
            if items:
                md.append(f"### {violation_type}")
                for item in items[:5]:  # 최대 5개만 표시
                    md.append(f"- {item}")
                if len(items) > 5:
                    md.append(f"- ... 외 {len(items) - 5}개")
                md.append("")

        return '\n'.join(md)


def main():
    """메인 함수"""
    import sys

    if len(sys.argv) < 2:
        print("Usage: python project_analyzer.py <project_root>")
        sys.exit(1)

    project_root = sys.argv[1]
    analyzer = ProjectAnalyzer(project_root)

    print("프로젝트 분석 중...")
    report = analyzer.analyze()

    print("\n" + analyzer.generate_report(report))


if __name__ == '__main__':
    main()
