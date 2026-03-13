#!/usr/bin/env python3
"""
Test script for PR auto-generation feature in skill-bug-reporter.

This script simulates the bug reporting workflow:
1. Create a feature branch
2. Generate a PR
3. Record the bug in skill_bugs.md with PR link
4. Verify the flow
"""

import subprocess
import json
from datetime import datetime
from pathlib import Path


def run_command(cmd, check=True):
    """Run a shell command and return output."""
    result = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
        check=False
    )
    if check and result.returncode != 0:
        print(f"❌ Command failed: {cmd}")
        print(f"   Error: {result.stderr}")
        return None
    return result.stdout.strip()


def test_pr_generation():
    """Test PR auto-generation workflow."""
    print("\n" + "=" * 60)
    print("Testing GitHub PR Auto-Generation Workflow")
    print("=" * 60 + "\n")

    # Setup
    repo_dir = Path(__file__).parent.parent.parent.parent.parent
    print(f"[*] Repository: {repo_dir}")

    # Test branch name
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    branch_name = f"fix/skill-bug-reporter-test-{timestamp}"

    print(f"[1] Creating test branch: {branch_name}")
    status = run_command(f"cd {repo_dir} && git branch -v", check=False)
    print(f"    Current branches:\n{status}\n")

    # Create feature branch
    create_branch = run_command(
        f"cd {repo_dir} && git checkout -b {branch_name}",
        check=False
    )
    if create_branch is None:
        print("[!] Branch creation might have failed (already exists?)")
    else:
        print(f"    ✅ Branch created: {branch_name}\n")

    # Simulate PR creation (dry-run)
    print("[2] Simulating PR creation (dry-run)...")
    pr_body = """## 오류 정보

- **스킬:** skill-bug-reporter
- **오류 유형:** test
- **오류 메시지:** Test PR auto-generation
- **발생 파일:** test_pr_auto_generation.py

## 예상 수정사항
Test PR for validating workflow integration.

자동 생성된 PR입니다. 준비 완료 후 main으로 merge합니다."""

    pr_cmd = (
        f'cd {repo_dir} && gh pr create '
        f'--title "[skill-bug-reporter] skill-bug-reporter: test 자동 수정" '
        f'--body "{pr_body}" '
        f'--draft 2>&1'
    )

    print(f"    PR Command (dry-run):")
    print(f"    {pr_cmd[:100]}...\n")

    # Attempt actual PR creation
    print("[3] Attempting actual PR creation...")
    pr_output = run_command(pr_cmd, check=False)

    if pr_output and "pull" in pr_output.lower():
        # Extract PR number from output
        # Output format: https://github.com/owner/repo/pull/NUMBER
        if "github.com" in pr_output and "pull" in pr_output:
            pr_parts = pr_output.split("/")
            pr_number = pr_parts[-1] if pr_parts else "UNKNOWN"
            print(f"    ✅ PR created successfully!")
            print(f"    PR URL: {pr_output}\n")

            # Record in skill_bugs.md
            print("[4] Recording in skill_bugs.md...")
            skill_bugs_path = Path.home() / ".claude" / "logs" / "skill_bugs.md"

            if not skill_bugs_path.exists():
                print(f"    Creating {skill_bugs_path}...")
                header = """# Skill Bug Log

| 날짜 | 스킬 | 유형 | 증상 요약 | 발생 파일 | 횟수 | PR 링크 | 상태 |
|------|------|------|----------|----------|------|--------|------|
"""
                skill_bugs_path.parent.mkdir(parents=True, exist_ok=True)
                skill_bugs_path.write_text(header)

            # Append test entry
            today = datetime.now().strftime("%Y-%m-%d")
            new_entry = f"| {today} | skill-bug-reporter | test | Test PR auto-generation | `test_pr_auto_generation.py` | 1회 | [PR #{pr_number}]({pr_output}) | PENDING |\n"

            with open(skill_bugs_path, "a") as f:
                f.write(new_entry)

            print(f"    ✅ Recorded in {skill_bugs_path}")
            print(f"    Entry: {new_entry.strip()}\n")

            return True
        else:
            print(f"    ⚠️  PR output format unexpected: {pr_output}\n")
            return False
    else:
        print(f"    ⚠️  PR creation failed or not authorized")
        print(f"    Output: {pr_output}\n")
        print("    Note: This is expected if GITHUB_TOKEN is not set.")
        print("    To enable PR creation, set: export GITHUB_TOKEN=<your-token>\n")

        # Show what would be created
        print("[*] Expected workflow (if authenticated):")
        print(f"    1. Branch: {branch_name}")
        print(f"    2. PR: [skill-bug-reporter] skill-bug-reporter: test 자동 수정")
        print(f"    3. skill_bugs.md updated with PR link")
        print(f"    4. On 2nd occurrence, Sonnet agent fixes + PR ready\n")

        return False


def cleanup_test_branch(branch_name):
    """Clean up test branch."""
    repo_dir = Path(__file__).parent.parent.parent.parent.parent
    print("[5] Cleanup (optional)...")
    print(f"    To delete test branch: git branch -D {branch_name}")
    print(f"    To delete remote: git push origin --delete {branch_name}\n")


if __name__ == "__main__":
    try:
        success = test_pr_generation()
        if success:
            print("\n✅ Test completed successfully!")
            print("   Next steps: Push branch + verify PR in GitHub\n")
        else:
            print("\n⚠️  Test completed with warnings.")
            print("   PR auto-creation requires GitHub authentication.\n")
    except Exception as e:
        print(f"\n❌ Test failed with error: {e}\n")
