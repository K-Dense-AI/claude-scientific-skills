"""
Asana Extended API Helper Module
MCP 도구에 없는 Asana REST API를 직접 호출하는 헬퍼.

사용법:
    from asana_api import AsanaAPI
    api = AsanaAPI()  # ASANA_PAT 환경변수에서 자동 로드
    api.add_task_to_project("TASK_GID", "PROJECT_GID")
"""

import os
import json
import time
import urllib.request
import urllib.error

BASE_URL = "https://app.asana.com/api/1.0"


class AsanaAPI:
    def __init__(self, pat=None):
        self.pat = pat or os.environ.get("ASANA_PAT")
        if not self.pat:
            raise ValueError(
                "ASANA_PAT 환경변수가 설정되지 않았습니다.\n"
                "1. https://app.asana.com/0/my-apps 접속\n"
                "2. Personal access token → New access token 생성\n"
                "3. 환경변수 설정: set ASANA_PAT=2/xxxx/xxxx:xxxx"
            )

    def _request(self, method, path, data=None):
        """Asana API 호출. 성공 시 응답 dict, 실패 시 예외 발생."""
        url = f"{BASE_URL}{path}"
        body = json.dumps({"data": data}).encode() if data else None

        req = urllib.request.Request(url, data=body, method=method, headers={
            "Authorization": f"Bearer {self.pat}",
            "Content-Type": "application/json",
        })

        try:
            resp = urllib.request.urlopen(req)
            content = resp.read().decode()
            return json.loads(content) if content.strip() else {}
        except urllib.error.HTTPError as e:
            error_body = e.read().decode()
            raise RuntimeError(f"Asana API 오류 {e.code}: {error_body[:500]}")

    # ── Task ↔ Project (multi-homing) ──

    def add_task_to_project(self, task_gid, project_gid, section_gid=None):
        """태스크를 프로젝트에 추가 (multi-homing)."""
        data = {"project": project_gid}
        if section_gid:
            data["section"] = section_gid
        return self._request("POST", f"/tasks/{task_gid}/addProject", data)

    def remove_task_from_project(self, task_gid, project_gid):
        """태스크를 프로젝트에서 제거."""
        return self._request("POST", f"/tasks/{task_gid}/removeProject",
                             {"project": project_gid})

    # ── Section 관리 ──

    def create_section(self, project_gid, name):
        """프로젝트에 섹션 생성. 생성된 섹션 dict 반환."""
        result = self._request("POST", f"/projects/{project_gid}/sections",
                               {"name": name})
        return result.get("data", {})

    def add_task_to_section(self, section_gid, task_gid):
        """태스크를 섹션으로 이동."""
        return self._request("POST", f"/sections/{section_gid}/addTask",
                             {"task": task_gid})

    def delete_section(self, section_gid):
        """빈 섹션 삭제."""
        return self._request("DELETE", f"/sections/{section_gid}")

    # ── Tag 관리 ──

    def create_tag(self, workspace_gid, name):
        """워크스페이스에 태그 생성."""
        result = self._request("POST", f"/workspaces/{workspace_gid}/tags",
                               {"name": name})
        return result.get("data", {})

    def add_tag_to_task(self, task_gid, tag_gid):
        """태스크에 태그 추가."""
        return self._request("POST", f"/tasks/{task_gid}/addTag",
                             {"tag": tag_gid})

    def remove_tag_from_task(self, task_gid, tag_gid):
        """태스크에서 태그 제거."""
        return self._request("POST", f"/tasks/{task_gid}/removeTag",
                             {"tag": tag_gid})

    # ── Portfolio 관리 ──

    def add_item_to_portfolio(self, portfolio_gid, item_gid):
        """포트폴리오에 항목(프로젝트) 추가."""
        return self._request("POST", f"/portfolios/{portfolio_gid}/addItem",
                             {"item": item_gid})

    def remove_item_from_portfolio(self, portfolio_gid, item_gid):
        """포트폴리오에서 항목 제거."""
        return self._request("POST", f"/portfolios/{portfolio_gid}/removeItem",
                             {"item": item_gid})

    # ── Project 복제 ──

    def duplicate_project(self, project_gid, new_name, include=None):
        """프로젝트 복제 (비동기). Job 정보 반환."""
        if include is None:
            include = ["members", "task_notes", "task_assignee",
                       "task_subtasks", "task_attachments", "task_dates",
                       "task_dependencies", "task_followers", "task_tags",
                       "task_projects"]
        result = self._request("POST", f"/projects/{project_gid}/duplicate",
                               {"name": new_name, "include": include})
        return result.get("data", {})

    # ── Batch 작업 ──

    def batch(self, actions):
        """여러 API 호출을 하나의 요청으로 묶음 (최대 ~10개)."""
        url = f"{BASE_URL}/batch"
        body = json.dumps({"data": {"actions": actions}}).encode()
        req = urllib.request.Request(url, data=body, method="POST", headers={
            "Authorization": f"Bearer {self.pat}",
            "Content-Type": "application/json",
        })
        resp = urllib.request.urlopen(req)
        return json.loads(resp.read().decode())

    def batch_add_to_project(self, task_gids, project_gid, section_gid=None):
        """여러 태스크를 한 프로젝트에 일괄 추가."""
        results = {"success": 0, "fail": 0, "details": []}

        # batch API는 ~10개 제한이므로 청크로 나눔
        for i in range(0, len(task_gids), 10):
            chunk = task_gids[i:i + 10]
            actions = []
            for gid in chunk:
                data = {"project": project_gid}
                if section_gid:
                    data["section"] = section_gid
                actions.append({
                    "method": "POST",
                    "relative_path": f"/tasks/{gid}/addProject",
                    "data": data,
                })

            try:
                resp = self.batch(actions)
                for item in resp.get("data", []):
                    status_code = item.get("status_code", 0)
                    if 200 <= status_code < 300:
                        results["success"] += 1
                    else:
                        results["fail"] += 1
                    results["details"].append(item)
            except Exception as e:
                # batch 실패 시 개별 호출로 fallback
                for gid in chunk:
                    try:
                        self.add_task_to_project(gid, project_gid, section_gid)
                        results["success"] += 1
                        time.sleep(0.3)
                    except Exception as e2:
                        results["fail"] += 1
                        results["details"].append({"task": gid, "error": str(e2)})

            if i + 10 < len(task_gids):
                time.sleep(0.5)

        return results


if __name__ == "__main__":
    # 간단 테스트
    api = AsanaAPI()
    print(f"PAT 로드 완료: ...{api.pat[-8:]}")
    print("사용 가능한 메서드:")
    for m in dir(api):
        if not m.startswith("_"):
            print(f"  api.{m}()")
