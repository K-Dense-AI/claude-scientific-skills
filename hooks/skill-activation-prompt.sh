#!/bin/bash
# Skill Activation Prompt Hook (UserPromptSubmit)
# stdin의 JSON을 Python 스크립트로 전달해 관련 스킬 제안 출력
cat | python '/c/Users/Jahyun/.claude/hooks/skill-activation-prompt.py'
