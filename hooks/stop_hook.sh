#!/bin/bash
# Stop hook: 응답 끝날 때 weekly + 현재 5h 블록 사용량 요약 저장

PENDING='/c/Users/Jahyun/.claude/pending_usage.txt'

TMPDIR='/c/Users/Jahyun/.claude'
ccusage weekly --offline --json 2>/dev/null > "$TMPDIR/claude_weekly.json"
ccusage blocks --recent --offline --json 2>/dev/null > "$TMPDIR/claude_blocks.json"

node << 'EOF' 2>/dev/null | tee "$PENDING"
const fs = require('fs');
const path = require('path');
const os = require('os');
const DIR = path.join(os.homedir(), '.claude');
try {
  const weekly = JSON.parse(fs.readFileSync(path.join(DIR, 'claude_weekly.json'), 'utf8'));
  const blocks = JSON.parse(fs.readFileSync(path.join(DIR, 'claude_blocks.json'), 'utf8'));

  const shortModel = name => {
    if (name.includes('opus'))   return 'opus';
    if (name.includes('sonnet')) return 'sonnet';
    if (name.includes('haiku'))  return 'haiku';
    return name.split('-')[1] || name;
  };
  const $$ = v => `$${Number(v).toFixed(0)}`;
  const pad = n => String(n).padStart(2, '0');

  // ── 5h 블록 ──
  const activeBlock = (blocks.blocks || []).find(b => b.isActive);
  let blockLine = '';
  if (activeBlock) {
    const TOKEN_LIMIT = 123657921;
    const curPct  = ((activeBlock.totalTokens / TOKEN_LIMIT) * 100).toFixed(0);
    const proj    = activeBlock.projection || {};
    const remainMin = proj.remainingMinutes || 0;
    const h = Math.floor(remainMin / 60);
    const m = pad(remainMin % 60);

    // 블록 초기화(종료) 시각
    const resetAt = new Date(Date.now() + remainMin * 60 * 1000);
    const resetStr = `${resetAt.getHours()}:${pad(resetAt.getMinutes())}`;

    // 예상 최종 토큰 % (현재 번레이트로 남은 시간 동안 추가 사용량 포함)
    const projTokens = proj.totalTokens || activeBlock.totalTokens;
    const projPct = ((projTokens / TOKEN_LIMIT) * 100).toFixed(0);
    const projCost = proj.totalCost ? $$(proj.totalCost) : '?';

    // 토큰 먼저 소진되는지 판단
    const tokenWarning = projTokens >= TOKEN_LIMIT
      ? `⚠ 토큰 ${projPct}% 초과 → 시간 전 차단 위험`
      : `여유 (예상최종 ${projPct}%)`;

    const blockPart = `블록${curPct}%·${$$(activeBlock.costUSD)}/${projCost}·리셋${resetStr}`;
    blockLine = projTokens >= TOKEN_LIMIT ? `[${blockPart}·!토큰초과]` : `[${blockPart}]`;
  }

  // ── 주간 ──
  const weekList = weekly.weekly || [];
  const currentWeek = weekList[weekList.length - 1];
  let weekLine = '';
  if (currentWeek) {
    const weekStart = new Date(currentWeek.week + 'T00:00:00');
    const weekEnd   = new Date(weekStart.getTime() + 7 * 24 * 3600 * 1000);
    const now = Date.now();

    // 경과율 (시간 기준)
    const weekPct = Math.min(100, (now - weekStart.getTime()) / (7 * 24 * 3600 * 1000) * 100);

    // 초기화까지 남은 시간
    const remainMs = weekEnd - now;
    const remainDays = Math.floor(remainMs / (24 * 3600 * 1000));
    const remainHrs  = Math.floor((remainMs % (24 * 3600 * 1000)) / (3600 * 1000));
    const resetStr = remainDays > 0 ? `D-${remainDays}` : `${remainHrs}h후`;

    // 역대 최대 대비 예상최종 %
    const pastWeeks = weekList.slice(0, -1);
    const maxPastCost = pastWeeks.length > 0 ? Math.max(...pastWeeks.map(w => w.totalCost)) : null;
    const projWeekRaw = weekPct > 0 ? currentWeek.totalCost / (weekPct / 100) : null;

    let overStr = '';
    if (maxPastCost && projWeekRaw) {
      const projPctNum = (projWeekRaw / maxPastCost) * 100;
      const projPct = projPctNum.toFixed(0);
      overStr = projPctNum >= 100 ? `·!예상${projPct}%` : `·예상${projPct}%`;
    }

    // 모델별
    const models = (currentWeek.modelBreakdowns || [])
      .filter(m => m.cost > 0.5)
      .sort((a, b) => b.cost - a.cost)
      .map(m => `${shortModel(m.modelName)}${$$(m.cost)}`)
      .join(' ');

    weekLine = `[주간${$$(currentWeek.totalCost)}·리셋${resetStr}${overStr}·${models}]`;
  }

  console.log([blockLine, weekLine].filter(Boolean).join(' '));
} catch(e) {
  console.log('사용량 파싱 오류: ' + e.message);
}
EOF
