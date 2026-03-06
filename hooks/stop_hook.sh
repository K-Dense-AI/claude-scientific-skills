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

const safeReadJSON = (filePath) => {
  try {
    const content = fs.readFileSync(filePath, 'utf8').trim();
    if (!content) return null;
    return JSON.parse(content);
  } catch(e) {
    return null;
  }
};

try {
  const weekly = safeReadJSON(path.join(DIR, 'claude_weekly.json'));
  const blocks = safeReadJSON(path.join(DIR, 'claude_blocks.json'));

  const shortModel = name => {
    if (name.includes('opus'))   return 'opus';
    if (name.includes('sonnet')) return 'sonnet';
    if (name.includes('haiku'))  return 'haiku';
    return name.split('-')[1] || name;
  };
  const pad = n => String(n).padStart(2, '0');

  // ── 5h 블록 ──
  const activeBlock = blocks ? (blocks.blocks || []).find(b => b.isActive) : null;
  let blockLine = '';
  if (activeBlock) {
    const TOKEN_LIMIT = 123657921;
    const curPct  = ((activeBlock.totalTokens / TOKEN_LIMIT) * 100).toFixed(0);
    const proj    = activeBlock.projection || {};
    const remainMin = proj.remainingMinutes || 0;
    const h = Math.floor(remainMin / 60);
    const m = pad(remainMin % 60);
    const resetStr = `${h}:${m}`;

    blockLine = `⏱ 5h | 토큰 ${curPct}% | 리셋까지 ${resetStr} 남음`;
  }

  // ── 주간 ──
  const weekList = weekly ? (weekly.weekly || []) : [];
  let weekLine = '';

  if (weekList.length === 0) {
    // 주간 리셋 직후 — 이번 주 데이터 없음
    weekLine = '📅 주간 | 방금 리셋됨 | 사용량 0%';
  } else {
    const currentWeek = weekList[weekList.length - 1];
    const weekStart = new Date(currentWeek.week + 'T00:00:00');
    const weekEnd   = new Date(weekStart.getTime() + 7 * 24 * 3600 * 1000);
    const now = Date.now();

    const weekPct = Math.min(100, (now - weekStart.getTime()) / (7 * 24 * 3600 * 1000) * 100);
    const remainMs = weekEnd - now;
    const remainDays = Math.floor(remainMs / (24 * 3600 * 1000));
    const remainHrs  = Math.floor((remainMs % (24 * 3600 * 1000)) / (3600 * 1000));
    const resetStr = remainDays > 0 ? `D-${remainDays}` : `${remainHrs}h후`;

    // 역대 최대 대비 예상최종 %
    const pastWeeks = weekList.slice(0, -1);
    const maxPastCost = pastWeeks.length > 0 ? Math.max(...pastWeeks.map(w => w.totalCost)) : null;
    const projWeekRaw = weekPct > 0 ? currentWeek.totalCost / (weekPct / 100) : null;
    let projStr = '';
    if (maxPastCost && projWeekRaw) {
      const projPct = ((projWeekRaw / maxPastCost) * 100).toFixed(0);
      projStr = ` | 예상 ${projPct}%`;
    }

    // 모델별 비율 (cost 기준 %)
    const breakdowns = (currentWeek.modelBreakdowns || []).filter(m => m.cost > 0.5);
    const totalCost = breakdowns.reduce((s, m) => s + m.cost, 0);
    const models = breakdowns
      .sort((a, b) => b.cost - a.cost)
      .map(m => `${shortModel(m.modelName)} ${totalCost > 0 ? ((m.cost / totalCost) * 100).toFixed(0) : 0}%`)
      .join('  ');

    weekLine = `📅 주간 | ${resetStr}${projStr}${models ? ' | ' + models : ''}`;
  }

  const lines = [blockLine, weekLine].filter(Boolean);
  console.log(lines.join('\n'));
} catch(e) {
  console.log('사용량 파싱 오류: ' + e.message);
}
EOF
