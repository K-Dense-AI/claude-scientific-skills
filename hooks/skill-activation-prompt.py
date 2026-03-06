#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Skill Activation Prompt Hook
UserPromptSubmit: analyze prompt and suggest relevant skills.
Based on diet103/claude-code-infrastructure-showcase pattern.
Output is ASCII-only for Windows cp949 compatibility.
"""
import io
import json
import re
import sys
from pathlib import Path


def load_rules():
    rules_path = Path.home() / '.claude' / 'skills' / 'skill-rules.json'
    if not rules_path.exists():
        return None
    with open(rules_path, encoding='utf-8') as f:
        return json.load(f)


def match_skill(prompt_lower, skill_cfg):
    triggers = skill_cfg.get('promptTriggers', {})

    for kw in triggers.get('keywords', []):
        if kw.lower() in prompt_lower:
            return 'keyword'

    for pat in triggers.get('intentPatterns', []):
        try:
            if re.search(pat, prompt_lower, re.IGNORECASE):
                return 'intent'
        except re.error:
            pass

    return None


def main():
    # Force UTF-8 output to avoid cp949 issues on Windows
    out = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

    try:
        raw = sys.stdin.buffer.read().decode('utf-8', errors='replace')
        data = json.loads(raw)
        prompt = data.get('prompt', '').lower()

        rules = load_rules()
        if not rules:
            sys.exit(0)

        skills = rules.get('skills', {})
        matched = {}

        for skill_name, cfg in skills.items():
            match_type = match_skill(prompt, cfg)
            if match_type:
                priority = cfg.get('priority', 'low')
                matched.setdefault(priority, []).append(skill_name)

        if not matched:
            sys.exit(0)

        lines = [
            '======================================',
            '[SKILL ACTIVATION]',
            '======================================',
            '',
        ]

        labels = {
            'critical': '[REQUIRED]',
            'high':     '[RECOMMENDED]',
            'medium':   '[SUGGESTED]',
            'low':      '[OPTIONAL]',
        }
        for priority in ('critical', 'high', 'medium', 'low'):
            if priority in matched:
                lines.append(labels[priority])
                for name in matched[priority]:
                    lines.append('  -> ' + name)
                lines.append('')

        lines.append('ACTION: Use Skill tool with the above skill(s) before responding')
        lines.append('======================================')

        out.write('\n'.join(lines) + '\n')
        out.flush()

    except Exception as e:
        sys.stderr.write('skill-activation-prompt error: ' + str(e) + '\n')

    sys.exit(0)


if __name__ == '__main__':
    main()
