import {
  Zap,
  FileText,
  TrendingUp,
  Clock,
  Plus,
  Sparkles,
  BarChart3,
} from "lucide-react";
import Link from "next/link";

const stats = [
  { label: "Active Campaigns", value: "3", icon: Zap, color: "text-purple-400", bg: "bg-purple-500/10" },
  { label: "Content Pieces", value: "47", icon: FileText, color: "text-blue-400", bg: "bg-blue-500/10" },
  { label: "Engagement Rate", value: "4.8%", icon: TrendingUp, color: "text-green-400", bg: "bg-green-500/10" },
  { label: "Scheduled Posts", value: "12", icon: Clock, color: "text-amber-400", bg: "bg-amber-500/10" },
];

const campaigns = [
  {
    name: "BrewCraft Summer Launch",
    status: "Active",
    platforms: ["Instagram", "TikTok"],
    posts: 18,
    engagement: "5.2%",
    statusColor: "bg-green-500",
  },
  {
    name: "Q1 Brand Awareness",
    status: "Active",
    platforms: ["Instagram", "Facebook", "LinkedIn"],
    posts: 24,
    engagement: "4.1%",
    statusColor: "bg-green-500",
  },
  {
    name: "Product Launch - Cold Brew",
    status: "Draft",
    platforms: ["Instagram", "TikTok"],
    posts: 0,
    engagement: "-",
    statusColor: "bg-amber-500",
  },
  {
    name: "Holiday Special 2025",
    status: "Completed",
    platforms: ["Instagram", "Facebook"],
    posts: 32,
    engagement: "6.3%",
    statusColor: "bg-gray-500",
  },
];

const quickActions = [
  {
    label: "New Campaign",
    description: "Create an AI-powered campaign from scratch",
    icon: Plus,
    href: "/campaign/new",
    color: "from-purple-600 to-purple-800",
  },
  {
    label: "Generate Content",
    description: "Use AI to create posts, reels, and stories",
    icon: Sparkles,
    href: "/content",
    color: "from-pink-600 to-pink-800",
  },
  {
    label: "View Analytics",
    description: "Track performance and engagement metrics",
    icon: BarChart3,
    href: "/analytics",
    color: "from-cyan-600 to-cyan-800",
  },
];

export default function DashboardPage() {
  return (
    <div className="p-8">
      {/* Header */}
      <div className="mb-8">
        <h1 className="text-3xl font-bold text-white">Campaign Dashboard</h1>
        <p className="text-gray-400 mt-1">
          Welcome back. Here&apos;s what&apos;s happening with your campaigns.
        </p>
      </div>

      {/* Stats Cards */}
      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-6 mb-8">
        {stats.map((stat) => {
          const Icon = stat.icon;
          return (
            <div
              key={stat.label}
              className="bg-gray-900 border border-gray-800 rounded-2xl p-6"
            >
              <div className="flex items-center justify-between mb-4">
                <div className={`${stat.bg} p-3 rounded-xl`}>
                  <Icon className={`w-5 h-5 ${stat.color}`} />
                </div>
              </div>
              <p className="text-3xl font-bold text-white">{stat.value}</p>
              <p className="text-sm text-gray-400 mt-1">{stat.label}</p>
            </div>
          );
        })}
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
        {/* Recent Campaigns */}
        <div className="lg:col-span-2 bg-gray-900 border border-gray-800 rounded-2xl p-6">
          <h2 className="text-lg font-semibold text-white mb-4">
            Recent Campaigns
          </h2>
          <div className="space-y-3">
            {campaigns.map((campaign) => (
              <div
                key={campaign.name}
                className="flex items-center justify-between p-4 bg-gray-800/50 rounded-xl hover:bg-gray-800 transition-colors"
              >
                <div className="flex items-center gap-4">
                  <div
                    className={`w-2.5 h-2.5 rounded-full ${campaign.statusColor}`}
                  />
                  <div>
                    <p className="text-sm font-medium text-white">
                      {campaign.name}
                    </p>
                    <p className="text-xs text-gray-400">
                      {campaign.platforms.join(" · ")}
                    </p>
                  </div>
                </div>
                <div className="flex items-center gap-6">
                  <div className="text-right">
                    <p className="text-sm text-white">{campaign.posts} posts</p>
                    <p className="text-xs text-gray-400">
                      {campaign.engagement} engagement
                    </p>
                  </div>
                  <span
                    className={`text-xs px-2.5 py-1 rounded-full font-medium ${
                      campaign.status === "Active"
                        ? "bg-green-500/10 text-green-400"
                        : campaign.status === "Draft"
                        ? "bg-amber-500/10 text-amber-400"
                        : "bg-gray-500/10 text-gray-400"
                    }`}
                  >
                    {campaign.status}
                  </span>
                </div>
              </div>
            ))}
          </div>
        </div>

        {/* Quick Actions */}
        <div className="space-y-4">
          <h2 className="text-lg font-semibold text-white">Quick Actions</h2>
          {quickActions.map((action) => {
            const Icon = action.icon;
            return (
              <Link
                key={action.label}
                href={action.href}
                className={`block p-5 rounded-2xl bg-gradient-to-br ${action.color} hover:opacity-90 transition-opacity`}
              >
                <Icon className="w-6 h-6 text-white mb-3" />
                <p className="font-semibold text-white">{action.label}</p>
                <p className="text-sm text-white/70 mt-1">
                  {action.description}
                </p>
              </Link>
            );
          })}
        </div>
      </div>
    </div>
  );
}
