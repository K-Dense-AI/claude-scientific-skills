#!/usr/bin/env python3
"""
PostToolUse hook: Agent 도구 완료 시 Notion에 로깅.
stdin으로 hook JSON 수신.
"""
import json, sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))


def main():
    try:
        hook_data = json.load(sys.stdin)
    except Exception:
        return

    tool_name = hook_data.get('tool_name', '')
    if tool_name != 'Agent':
        return

    tool_input = hook_data.get('tool_input', {})
    tool_response = hook_data.get('tool_response', {})

    # 에이전트 정보 추출
    agent_name = tool_input.get('name', 'unknown-agent')
    agent_desc = tool_input.get('description', '')
    team_name  = tool_input.get('team_name', '')
    prompt     = tool_input.get('prompt', '')[:300]

    # 응답에서 결과 및 토큰 추출
    response_content = ''
    tokens_in = None
    tokens_out = None
    cost = None
    elapsed = None
    if isinstance(tool_response, dict):
        response_content = str(tool_response.get('output', '') or tool_response.get('result', ''))[:500]
        usage = tool_response.get('usage', {}) or {}
        if isinstance(usage, dict) and usage:
            tokens_in = int(usage.get('input_tokens', 0) or 0) or None
            tokens_out = int(usage.get('output_tokens', 0) or 0) or None
        cost_val = tool_response.get('cost_usd') or tool_response.get('total_cost_usd')
        cost = float(cost_val) if cost_val is not None else None
        dur_ms = tool_response.get('duration_ms')
        elapsed = round(dur_ms / 1000, 1) if dur_ms else None
    elif isinstance(tool_response, str):
        response_content = tool_response[:500]

    # project_id = team_name (없으면 'sdk-direct')
    project_id = team_name or 'sdk-direct'

    try:
        import notion_logger as nl
        nl.log_event(
            event='agent_done',
            name=agent_name,
            agent_id=f'sdk-{agent_name}',
            project_id=project_id,
            team_type='general',
            role='solo',
            status='done',
            task=agent_desc or prompt,
            result=response_content,
            tokens_in=tokens_in,
            tokens_out=tokens_out,
            cost=cost,
            elapsed=elapsed,
        )
        print(f'[notion-hook] logged: {agent_name} / {project_id}', flush=True)
    except Exception as e:
        print(f'[notion-hook] error: {e}', flush=True)


if __name__ == '__main__':
    main()
