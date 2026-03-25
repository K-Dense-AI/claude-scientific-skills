"use client";

import { useState } from "react";
import { Eye, TrendingUp, Users, MousePointerClick, Lightbulb } from "lucide-react";
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  Area,
  AreaChart,
} from "recharts";

const dateRanges = ["7d", "14d", "30d", "90d"];

const engagementData7d = [
  { date: "Mon", engagement: 4.2, reach: 3200 },
  { date: "Tue", engagement: 4.8, reach: 4100 },
  { date: "Wed", engagement: 3.9, reach: 2800 },
  { date: "Thu", engagement: 5.1, reach: 5200 },
  { date: "Fri", engagement: 5.6, reach: 6100 },
  { date: "Sat", engagement: 4.4, reach: 3800 },
  { date: "Sun", engagement: 4.9, reach: 4500 },
];

const engagementData30d = [
  { date: "Week 1", engagement: 4.1, reach: 22000 },
  { date: "Week 2", engagement: 4.5, reach: 28000 },
  { date: "Week 3", engagement: 4.8, reach: 31000 },
  { date: "Week 4", engagement: 5.2, reach: 35000 },
];

const platformData = [
  { platform: "Instagram", followers: 12400, engagement: 5.2, reach: 45000, color: "#ec4899" },
  { platform: "Facebook", followers: 8200, engagement: 3.8, reach: 28000, color: "#3b82f6" },
  { platform: "LinkedIn", followers: 4500, engagement: 4.1, reach: 15000, color: "#06b6d4" },
  { platform: "TikTok", followers: 18900, engagement: 7.3, reach: 82000, color: "#ef4444" },
];

const topPosts = [
  {
    caption: "POV: You just discovered the perfect pour-over ratio...",
    platform: "TikTok",
    reach: "82K",
    engagement: "7.3%",
    likes: "3.1K",
  },
  {
    caption: "Discover our new single-origin beans from Ethiopia...",
    platform: "Instagram",
    reach: "45K",
    engagement: "5.2%",
    likes: "892",
  },
  {
    caption: "How we reduced our carbon footprint by 40%...",
    platform: "LinkedIn",
    reach: "15K",
    engagement: "4.1%",
    likes: "423",
  },
];

const recommendations = [
  {
    title: "Post More Reels",
    description:
      "Your Reels get 3.2x more engagement than static posts. Aim for 4-5 Reels per week.",
  },
  {
    title: "Optimal Posting Times",
    description:
      "Your audience is most active between 9-11 AM and 7-9 PM. Schedule key posts during these windows.",
  },
  {
    title: "Leverage TikTok Growth",
    description:
      "TikTok is your fastest-growing platform (+23% this month). Increase content frequency from 3 to 5 posts/week.",
  },
];

const kpis = [
  { label: "Total Reach", value: "170K", change: "+12%", icon: Eye, color: "text-purple-400", bg: "bg-purple-500/10" },
  { label: "Engagement Rate", value: "4.8%", change: "+0.6%", icon: TrendingUp, color: "text-green-400", bg: "bg-green-500/10" },
  { label: "Followers Growth", value: "+2.4K", change: "+8%", icon: Users, color: "text-blue-400", bg: "bg-blue-500/10" },
  { label: "Click-Through Rate", value: "2.3%", change: "+0.4%", icon: MousePointerClick, color: "text-amber-400", bg: "bg-amber-500/10" },
];

export default function AnalyticsPage() {
  const [range, setRange] = useState("7d");
  const chartData = range === "30d" || range === "90d" ? engagementData30d : engagementData7d;

  return (
    <div className="p-8">
      <div className="flex items-center justify-between mb-8">
        <div>
          <h1 className="text-3xl font-bold text-white">Analytics</h1>
          <p className="text-gray-400 mt-1">
            Track campaign performance and audience growth.
          </p>
        </div>
        <div className="flex gap-1 bg-gray-900 rounded-xl p-1 border border-gray-800">
          {dateRanges.map((r) => (
            <button
              key={r}
              onClick={() => setRange(r)}
              className={`px-4 py-2 rounded-lg text-sm font-medium transition-colors ${
                range === r
                  ? "bg-purple-600 text-white"
                  : "text-gray-400 hover:text-white"
              }`}
            >
              {r}
            </button>
          ))}
        </div>
      </div>

      {/* KPI Cards */}
      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-6 mb-8">
        {kpis.map((kpi) => {
          const Icon = kpi.icon;
          return (
            <div
              key={kpi.label}
              className="bg-gray-900 border border-gray-800 rounded-2xl p-6"
            >
              <div className="flex items-center justify-between mb-4">
                <div className={`${kpi.bg} p-3 rounded-xl`}>
                  <Icon className={`w-5 h-5 ${kpi.color}`} />
                </div>
                <span className="text-xs font-medium text-green-400 bg-green-500/10 px-2 py-1 rounded-full">
                  {kpi.change}
                </span>
              </div>
              <p className="text-3xl font-bold text-white">{kpi.value}</p>
              <p className="text-sm text-gray-400 mt-1">{kpi.label}</p>
            </div>
          );
        })}
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 mb-8">
        {/* Engagement Line Chart */}
        <div className="bg-gray-900 border border-gray-800 rounded-2xl p-6">
          <h2 className="text-lg font-semibold text-white mb-4">
            Engagement Over Time
          </h2>
          <ResponsiveContainer width="100%" height={280}>
            <AreaChart data={chartData}>
              <defs>
                <linearGradient id="engGrad" x1="0" y1="0" x2="0" y2="1">
                  <stop offset="5%" stopColor="#a855f7" stopOpacity={0.3} />
                  <stop offset="95%" stopColor="#a855f7" stopOpacity={0} />
                </linearGradient>
              </defs>
              <CartesianGrid strokeDasharray="3 3" stroke="#1f2937" />
              <XAxis dataKey="date" stroke="#6b7280" fontSize={12} />
              <YAxis stroke="#6b7280" fontSize={12} />
              <Tooltip
                contentStyle={{
                  backgroundColor: "#1f2937",
                  border: "1px solid #374151",
                  borderRadius: "8px",
                  color: "#fff",
                  fontSize: "12px",
                }}
              />
              <Area
                type="monotone"
                dataKey="engagement"
                stroke="#a855f7"
                strokeWidth={2}
                fill="url(#engGrad)"
              />
            </AreaChart>
          </ResponsiveContainer>
        </div>

        {/* Platform Bar Chart */}
        <div className="bg-gray-900 border border-gray-800 rounded-2xl p-6">
          <h2 className="text-lg font-semibold text-white mb-4">
            Performance by Platform
          </h2>
          <ResponsiveContainer width="100%" height={280}>
            <BarChart data={platformData}>
              <CartesianGrid strokeDasharray="3 3" stroke="#1f2937" />
              <XAxis dataKey="platform" stroke="#6b7280" fontSize={12} />
              <YAxis stroke="#6b7280" fontSize={12} />
              <Tooltip
                contentStyle={{
                  backgroundColor: "#1f2937",
                  border: "1px solid #374151",
                  borderRadius: "8px",
                  color: "#fff",
                  fontSize: "12px",
                }}
              />
              <Bar dataKey="reach" radius={[6, 6, 0, 0]} fill="#a855f7" />
            </BarChart>
          </ResponsiveContainer>
        </div>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Top Performing Posts */}
        <div className="bg-gray-900 border border-gray-800 rounded-2xl p-6">
          <h2 className="text-lg font-semibold text-white mb-4">
            Top Performing Posts
          </h2>
          <div className="space-y-3">
            {topPosts.map((post, i) => (
              <div
                key={i}
                className="flex items-start gap-4 p-4 bg-gray-800/50 rounded-xl"
              >
                <div className="flex items-center justify-center w-8 h-8 rounded-full bg-purple-500/10 text-purple-400 text-sm font-bold shrink-0">
                  {i + 1}
                </div>
                <div className="flex-1 min-w-0">
                  <p className="text-sm text-white truncate">{post.caption}</p>
                  <div className="flex items-center gap-3 mt-1.5 text-xs text-gray-500">
                    <span>{post.platform}</span>
                    <span>Reach: {post.reach}</span>
                    <span>Eng: {post.engagement}</span>
                    <span>Likes: {post.likes}</span>
                  </div>
                </div>
              </div>
            ))}
          </div>
        </div>

        {/* Growth Recommendations */}
        <div className="bg-gray-900 border border-gray-800 rounded-2xl p-6">
          <h2 className="text-lg font-semibold text-white mb-4">
            Growth Recommendations
          </h2>
          <div className="space-y-3">
            {recommendations.map((rec, i) => (
              <div
                key={i}
                className="flex items-start gap-3 p-4 bg-gray-800/50 rounded-xl"
              >
                <div className="p-2 bg-amber-500/10 rounded-lg shrink-0">
                  <Lightbulb className="w-4 h-4 text-amber-400" />
                </div>
                <div>
                  <p className="text-sm font-medium text-white">{rec.title}</p>
                  <p className="text-xs text-gray-400 mt-1">
                    {rec.description}
                  </p>
                </div>
              </div>
            ))}
          </div>
        </div>
      </div>
    </div>
  );
}
