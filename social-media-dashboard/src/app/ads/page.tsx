"use client";

import { useState } from "react";
import {
  Megaphone,
  Target,
  DollarSign,
  TrendingUp,
  Eye,
  MousePointerClick,
  ExternalLink,
  Zap,
  Globe,
  Send,
} from "lucide-react";

const adPlatforms = [
  {
    id: "meta",
    name: "Meta Ads",
    subtitle: "Facebook & Instagram",
    color: "bg-blue-500",
    textColor: "text-blue-400",
    borderColor: "border-blue-500/30",
    tools: 20,
    capabilities: ["Image Campaigns", "Carousel Ads", "Audience Targeting", "Retargeting", "Lookalike Audiences"],
    status: "connected" as const,
  },
  {
    id: "google",
    name: "Google Ads",
    subtitle: "Search & Display",
    color: "bg-yellow-500",
    textColor: "text-yellow-400",
    borderColor: "border-yellow-500/30",
    tools: 39,
    capabilities: ["Keyword Research", "Search Campaigns", "Performance Max", "Ad Extensions", "Smart Bidding"],
    status: "connected" as const,
  },
  {
    id: "linkedin",
    name: "LinkedIn Ads",
    subtitle: "B2B Advertising",
    color: "bg-cyan-500",
    textColor: "text-cyan-400",
    borderColor: "border-cyan-500/30",
    tools: 28,
    capabilities: ["Sponsored Content", "Lead Gen Forms", "Account Targeting", "InMail Campaigns", "Analytics"],
    status: "connected" as const,
  },
  {
    id: "tiktok",
    name: "TikTok Ads",
    subtitle: "In-Feed Campaigns",
    color: "bg-red-500",
    textColor: "text-red-400",
    borderColor: "border-red-500/30",
    tools: 4,
    capabilities: ["In-Feed Video Ads", "Spark Ads", "TopView Ads", "Asset Validation"],
    status: "connected" as const,
  },
];

const postizPlatforms = [
  { name: "Instagram", icon: "📸", posts: 12 },
  { name: "Facebook", icon: "📘", posts: 8 },
  { name: "LinkedIn", icon: "💼", posts: 15 },
  { name: "TikTok", icon: "🎵", posts: 6 },
  { name: "Twitter/X", icon: "🐦", posts: 10 },
  { name: "YouTube", icon: "🎬", posts: 3 },
  { name: "Pinterest", icon: "📌", posts: 4 },
  { name: "Reddit", icon: "🤖", posts: 2 },
];

const activeCampaigns = [
  {
    name: "Spring Collection Launch",
    platform: "Meta Ads",
    budget: "$2,500",
    spent: "$1,847",
    impressions: "125.4K",
    clicks: "3,289",
    ctr: "2.62%",
    conversions: 47,
    status: "active",
  },
  {
    name: "Brand Awareness Q1",
    platform: "Google Ads",
    budget: "$5,000",
    spent: "$3,210",
    impressions: "890.2K",
    clicks: "12,450",
    ctr: "1.40%",
    conversions: 156,
    status: "active",
  },
  {
    name: "B2B Lead Gen",
    platform: "LinkedIn Ads",
    budget: "$3,000",
    spent: "$2,100",
    impressions: "45.8K",
    clicks: "892",
    ctr: "1.95%",
    conversions: 23,
    status: "active",
  },
];

export default function AdsPage() {
  const [activeTab, setActiveTab] = useState<"platforms" | "campaigns" | "scheduler">("platforms");

  return (
    <div className="p-8 space-y-8">
      <div>
        <h1 className="text-3xl font-bold text-white flex items-center gap-3">
          <Megaphone className="w-8 h-8 text-purple-400" />
          Ads & Publishing
        </h1>
        <p className="text-gray-400 mt-1">
          Manage paid campaigns via Adspirer (91 tools) and organic publishing via Postiz (28+ platforms)
        </p>
      </div>

      {/* Tabs */}
      <div className="flex gap-2 border-b border-gray-800 pb-0">
        {[
          { id: "platforms" as const, label: "Ad Platforms", icon: Target },
          { id: "campaigns" as const, label: "Active Campaigns", icon: TrendingUp },
          { id: "scheduler" as const, label: "Organic Publishing", icon: Send },
        ].map((tab) => (
          <button
            key={tab.id}
            onClick={() => setActiveTab(tab.id)}
            className={`flex items-center gap-2 px-5 py-3 text-sm font-medium rounded-t-xl transition-colors ${
              activeTab === tab.id
                ? "bg-gray-900 text-white border border-gray-800 border-b-gray-900 -mb-px"
                : "text-gray-500 hover:text-gray-300"
            }`}
          >
            <tab.icon className="w-4 h-4" />
            {tab.label}
          </button>
        ))}
      </div>

      {/* Ad Platforms Tab */}
      {activeTab === "platforms" && (
        <div className="space-y-6">
          {/* Adspirer Integration */}
          <div className="flex items-center justify-between">
            <div>
              <h2 className="text-lg font-semibold text-white">Adspirer — 91 Ad Management Tools</h2>
              <p className="text-sm text-gray-500">Connected via MCP Server • OAuth 2.1 authenticated</p>
            </div>
            <div className="flex items-center gap-2 px-3 py-1.5 bg-green-500/10 text-green-400 rounded-full text-xs font-medium">
              <div className="w-2 h-2 rounded-full bg-green-400 animate-pulse" />
              Connected
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            {adPlatforms.map((platform) => (
              <div
                key={platform.id}
                className={`bg-gray-900 border ${platform.borderColor} rounded-2xl p-5 space-y-4`}
              >
                <div className="flex items-center justify-between">
                  <div className="flex items-center gap-3">
                    <div className={`w-10 h-10 ${platform.color} rounded-xl flex items-center justify-center`}>
                      <Globe className="w-5 h-5 text-white" />
                    </div>
                    <div>
                      <h3 className="text-white font-semibold">{platform.name}</h3>
                      <p className="text-xs text-gray-500">{platform.subtitle}</p>
                    </div>
                  </div>
                  <span className="text-xs bg-gray-800 text-gray-400 px-2.5 py-1 rounded-full">
                    {platform.tools} tools
                  </span>
                </div>

                <div className="flex flex-wrap gap-1.5">
                  {platform.capabilities.map((cap) => (
                    <span
                      key={cap}
                      className="text-xs px-2 py-1 bg-gray-800 text-gray-400 rounded-lg"
                    >
                      {cap}
                    </span>
                  ))}
                </div>

                <div className="flex gap-2">
                  <button className={`flex-1 py-2 px-3 ${platform.color} hover:opacity-90 text-white text-sm font-medium rounded-xl transition-opacity flex items-center justify-center gap-1.5`}>
                    <Zap className="w-3.5 h-3.5" />
                    Create Campaign
                  </button>
                  <button className="py-2 px-3 bg-gray-800 hover:bg-gray-700 text-gray-300 text-sm rounded-xl transition-colors flex items-center gap-1.5">
                    <ExternalLink className="w-3.5 h-3.5" />
                    Dashboard
                  </button>
                </div>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Active Campaigns Tab */}
      {activeTab === "campaigns" && (
        <div className="space-y-4">
          <div className="grid grid-cols-4 gap-4 mb-6">
            {[
              { label: "Total Budget", value: "$10,500", icon: DollarSign, color: "text-green-400" },
              { label: "Total Spent", value: "$7,157", icon: TrendingUp, color: "text-purple-400" },
              { label: "Total Impressions", value: "1.06M", icon: Eye, color: "text-blue-400" },
              { label: "Total Clicks", value: "16,631", icon: MousePointerClick, color: "text-yellow-400" },
            ].map((stat) => (
              <div key={stat.label} className="bg-gray-900 border border-gray-800 rounded-2xl p-4">
                <div className="flex items-center gap-2 mb-2">
                  <stat.icon className={`w-4 h-4 ${stat.color}`} />
                  <span className="text-xs text-gray-500">{stat.label}</span>
                </div>
                <p className="text-2xl font-bold text-white">{stat.value}</p>
              </div>
            ))}
          </div>

          <div className="bg-gray-900 border border-gray-800 rounded-2xl overflow-hidden">
            <table className="w-full">
              <thead>
                <tr className="border-b border-gray-800">
                  {["Campaign", "Platform", "Budget", "Spent", "Impressions", "Clicks", "CTR", "Conv.", "Status"].map((h) => (
                    <th key={h} className="text-left text-xs font-medium text-gray-500 px-4 py-3 uppercase tracking-wider">{h}</th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {activeCampaigns.map((c) => (
                  <tr key={c.name} className="border-b border-gray-800/50 hover:bg-gray-800/30 transition-colors">
                    <td className="px-4 py-3 text-sm text-white font-medium">{c.name}</td>
                    <td className="px-4 py-3 text-sm text-gray-400">{c.platform}</td>
                    <td className="px-4 py-3 text-sm text-gray-400">{c.budget}</td>
                    <td className="px-4 py-3 text-sm text-white">{c.spent}</td>
                    <td className="px-4 py-3 text-sm text-gray-400">{c.impressions}</td>
                    <td className="px-4 py-3 text-sm text-gray-400">{c.clicks}</td>
                    <td className="px-4 py-3 text-sm text-purple-400 font-medium">{c.ctr}</td>
                    <td className="px-4 py-3 text-sm text-green-400 font-medium">{c.conversions}</td>
                    <td className="px-4 py-3">
                      <span className="px-2 py-1 bg-green-500/10 text-green-400 rounded-full text-xs font-medium">Active</span>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {/* Organic Publishing Tab (Postiz) */}
      {activeTab === "scheduler" && (
        <div className="space-y-6">
          <div className="flex items-center justify-between">
            <div>
              <h2 className="text-lg font-semibold text-white">Postiz — Organic Publishing</h2>
              <p className="text-sm text-gray-500">Schedule and publish across 28+ social platforms</p>
            </div>
            <div className="flex items-center gap-2 px-3 py-1.5 bg-green-500/10 text-green-400 rounded-full text-xs font-medium">
              <div className="w-2 h-2 rounded-full bg-green-400 animate-pulse" />
              Connected
            </div>
          </div>

          <div className="grid grid-cols-2 md:grid-cols-4 gap-3">
            {postizPlatforms.map((p) => (
              <div key={p.name} className="bg-gray-900 border border-gray-800 rounded-xl p-4 flex items-center gap-3">
                <span className="text-2xl">{p.icon}</span>
                <div>
                  <p className="text-sm text-white font-medium">{p.name}</p>
                  <p className="text-xs text-gray-500">{p.posts} scheduled</p>
                </div>
              </div>
            ))}
          </div>

          <div className="bg-gray-900 border border-gray-800 rounded-2xl p-6 space-y-4">
            <h3 className="text-white font-semibold">Quick Publish</h3>
            <p className="text-sm text-gray-500">
              Create and schedule posts across all connected platforms. Supports images, videos, carousels, and threads.
            </p>
            <div className="grid grid-cols-1 md:grid-cols-3 gap-3">
              <button className="py-3 px-4 bg-purple-600 hover:bg-purple-700 text-white text-sm font-medium rounded-xl transition-colors flex items-center justify-center gap-2">
                <Send className="w-4 h-4" />
                New Post
              </button>
              <button className="py-3 px-4 bg-gray-800 hover:bg-gray-700 text-white text-sm font-medium rounded-xl transition-colors border border-gray-700 flex items-center justify-center gap-2">
                <Globe className="w-4 h-4" />
                Bulk Schedule
              </button>
              <button className="py-3 px-4 bg-gray-800 hover:bg-gray-700 text-white text-sm font-medium rounded-xl transition-colors border border-gray-700 flex items-center justify-center gap-2">
                <TrendingUp className="w-4 h-4" />
                View Analytics
              </button>
            </div>
          </div>

          <div className="bg-gray-900 border border-gray-800 rounded-2xl p-6 space-y-3">
            <h3 className="text-white font-semibold">Postiz Capabilities</h3>
            <div className="grid grid-cols-2 md:grid-cols-3 gap-2">
              {[
                "Multi-platform posting",
                "Media upload (PNG/JPG/MP4)",
                "Thread support (X/Twitter)",
                "Carousel posts (LinkedIn/IG)",
                "Story publishing",
                "Community posting (Reddit)",
                "Analytics per post",
                "Batch scheduling",
                "Platform-specific settings",
                "Draft management",
                "Post deletion",
                "Integration discovery",
              ].map((cap) => (
                <div key={cap} className="flex items-center gap-2 text-sm text-gray-400">
                  <div className="w-1.5 h-1.5 rounded-full bg-purple-500" />
                  {cap}
                </div>
              ))}
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
