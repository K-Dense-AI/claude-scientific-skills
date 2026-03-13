"""GitHub 최근 커밋 수집 - 주간 브리핑용"""
import json
import sys
import io
import urllib.request
import urllib.error
from datetime import datetime, timedelta, timezone
from pathlib import Path

# Windows cp949 환경에서 한글 출력 시 깨짐 방지
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

TEMP_DIR = Path(r"C:\Users\Jahyun\.claude\briefing_temp")
OUTPUT_FILE = TEMP_DIR / "briefing_github.json"

secrets = json.loads(Path("~/.claude/secrets.json").expanduser().read_text())
TOKEN = secrets["GITHUB_PAT"]
OWNER = "jahyunlee00299"
REPOS = ["UDH_Clustering", "Kinetic-modeling-and-optimization", "biosteam-tagatose"]

def gh_get(url):
    req = urllib.request.Request(url, headers={
        "Authorization": f"Bearer {TOKEN}",
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28"
    })
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read())

def collect():
    since = (datetime.now(timezone.utc) - timedelta(days=7)).strftime("%Y-%m-%dT%H:%M:%SZ")
    results = {"repos": [], "collected_at": str(datetime.now(timezone.utc).date())}

    for repo in REPOS:
        try:
            commits_url = f"https://api.github.com/repos/{OWNER}/{repo}/commits?since={since}&per_page=10"
            commits = gh_get(commits_url)
            repo_info = gh_get(f"https://api.github.com/repos/{OWNER}/{repo}")

            results["repos"].append({
                "name": repo,
                "url": repo_info.get("html_url"),
                "last_push": repo_info.get("pushed_at", "")[:10],
                "recent_commits": [
                    {
                        "sha": c["sha"][:7],
                        "message": c["commit"]["message"].split("\n")[0][:80],
                        "date": c["commit"]["author"]["date"][:10],
                        "author": c["commit"]["author"]["name"]
                    }
                    for c in commits[:5]
                ],
                "commit_count_7d": len(commits)
            })
        except Exception as e:
            results["repos"].append({"name": repo, "error": str(e)})

    TEMP_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_FILE.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    active = sum(1 for r in results["repos"] if r.get("commit_count_7d", 0) > 0)
    print(f"[GitHub] done: {len(REPOS)} repos, active={active}")
    print(f"saved: {OUTPUT_FILE}")

if __name__ == "__main__":
    collect()
