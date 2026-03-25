"use client";

import { useState } from "react";
import { Plus, Eye, Heart, MessageCircle } from "lucide-react";

type ContentPiece = {
  id: number;
  platform: string;
  type: string;
  caption: string;
  status: "Draft" | "Ready" | "Published";
  gradient: string;
  engagement?: { views: string; likes: string; comments: string };
};

const allContent: ContentPiece[] = [
  {
    id: 1,
    platform: "Instagram",
    type: "Carousel",
    caption:
      "Discover our new single-origin beans from Ethiopia. Rich, complex, and unforgettable flavors that tell a story of generations of farming.",
    status: "Published",
    gradient: "from-purple-600 to-pink-600",
    engagement: { views: "12.4K", likes: "892", comments: "67" },
  },
  {
    id: 2,
    platform: "TikTok",
    type: "Reel",
    caption:
      "POV: You just discovered the perfect pour-over ratio and your morning will never be the same again.",
    status: "Published",
    gradient: "from-red-600 to-orange-500",
    engagement: { views: "45.2K", likes: "3.1K", comments: "234" },
  },
  {
    id: 3,
    platform: "LinkedIn",
    type: "Post",
    caption:
      "How we reduced our carbon footprint by 40% this quarter while scaling production across 3 new markets.",
    status: "Published",
    gradient: "from-cyan-600 to-blue-600",
    engagement: { views: "8.7K", likes: "423", comments: "56" },
  },
  {
    id: 4,
    platform: "Instagram",
    type: "Reel",
    caption:
      "3 cold brew recipes you absolutely need to try this summer. The last one will surprise you!",
    status: "Ready",
    gradient: "from-pink-500 to-rose-600",
  },
  {
    id: 5,
    platform: "Facebook",
    type: "Post",
    caption:
      "Join us this Saturday for our free coffee cupping event at the downtown cafe. Limited spots available!",
    status: "Ready",
    gradient: "from-blue-500 to-indigo-600",
  },
  {
    id: 6,
    platform: "Instagram",
    type: "Story",
    caption:
      "Quick poll: Light roast or dark roast? Drop your vote below and we'll reveal the winner tomorrow!",
    status: "Ready",
    gradient: "from-violet-500 to-purple-700",
  },
  {
    id: 7,
    platform: "TikTok",
    type: "Reel",
    caption:
      "Stop doing pour-over coffee like this. Here's what the pros actually do differently.",
    status: "Draft",
    gradient: "from-red-500 to-pink-600",
  },
  {
    id: 8,
    platform: "Instagram",
    type: "Carousel",
    caption:
      "From bean to cup: the complete journey of your morning coffee. A story of passion, craft, and community.",
    status: "Draft",
    gradient: "from-amber-500 to-orange-600",
  },
  {
    id: 9,
    platform: "LinkedIn",
    type: "Post",
    caption:
      "Excited to announce our partnership with Fair Trade certified cooperatives across Central America.",
    status: "Draft",
    gradient: "from-teal-500 to-cyan-600",
  },
];

const tabs = ["All", "Posts", "Reels", "Carousels", "Stories"];

const platformBadgeColors: Record<string, string> = {
  Instagram: "bg-pink-500/20 text-pink-400",
  Facebook: "bg-blue-500/20 text-blue-400",
  LinkedIn: "bg-cyan-500/20 text-cyan-400",
  TikTok: "bg-red-500/20 text-red-400",
};

const statusColors: Record<string, string> = {
  Draft: "bg-gray-500/20 text-gray-400",
  Ready: "bg-amber-500/20 text-amber-400",
  Published: "bg-green-500/20 text-green-400",
};

export default function ContentPage() {
  const [activeTab, setActiveTab] = useState("All");

  const filtered =
    activeTab === "All"
      ? allContent
      : allContent.filter((c) => {
          if (activeTab === "Posts") return c.type === "Post";
          if (activeTab === "Reels") return c.type === "Reel";
          if (activeTab === "Carousels") return c.type === "Carousel";
          if (activeTab === "Stories") return c.type === "Story";
          return true;
        });

  return (
    <div className="p-8">
      <div className="flex items-center justify-between mb-8">
        <div>
          <h1 className="text-3xl font-bold text-white">Content Library</h1>
          <p className="text-gray-400 mt-1">
            Browse, edit, and manage all your content pieces.
          </p>
        </div>
        <button className="flex items-center gap-2 px-5 py-2.5 bg-purple-600 hover:bg-purple-700 text-white rounded-xl text-sm font-medium transition-colors">
          <Plus className="w-4 h-4" />
          Generate New
        </button>
      </div>

      {/* Filter Tabs */}
      <div className="flex gap-1 mb-6 bg-gray-900 rounded-xl p-1 w-fit border border-gray-800">
        {tabs.map((tab) => (
          <button
            key={tab}
            onClick={() => setActiveTab(tab)}
            className={`px-4 py-2 rounded-lg text-sm font-medium transition-colors ${
              activeTab === tab
                ? "bg-purple-600 text-white"
                : "text-gray-400 hover:text-white"
            }`}
          >
            {tab}
          </button>
        ))}
      </div>

      {/* Content Grid */}
      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
        {filtered.map((piece) => (
          <div
            key={piece.id}
            className="bg-gray-900 border border-gray-800 rounded-2xl overflow-hidden hover:border-gray-700 transition-colors"
          >
            {/* Visual Placeholder */}
            <div
              className={`h-48 bg-gradient-to-br ${piece.gradient} flex items-center justify-center`}
            >
              <span className="text-white/30 text-sm font-medium uppercase tracking-wider">
                {piece.type} Preview
              </span>
            </div>

            <div className="p-5 space-y-3">
              {/* Badges */}
              <div className="flex items-center gap-2">
                <span
                  className={`text-xs px-2.5 py-1 rounded-full font-medium ${platformBadgeColors[piece.platform]}`}
                >
                  {piece.platform}
                </span>
                <span className="text-xs px-2.5 py-1 rounded-full bg-gray-800 text-gray-300">
                  {piece.type}
                </span>
                <span
                  className={`text-xs px-2.5 py-1 rounded-full font-medium ml-auto ${statusColors[piece.status]}`}
                >
                  {piece.status}
                </span>
              </div>

              {/* Caption */}
              <p className="text-sm text-gray-300 leading-relaxed line-clamp-3">
                {piece.caption.slice(0, 100)}
                {piece.caption.length > 100 ? "..." : ""}
              </p>

              {/* Engagement */}
              {piece.engagement && (
                <div className="flex items-center gap-4 pt-2 border-t border-gray-800">
                  <div className="flex items-center gap-1 text-xs text-gray-500">
                    <Eye className="w-3.5 h-3.5" />
                    {piece.engagement.views}
                  </div>
                  <div className="flex items-center gap-1 text-xs text-gray-500">
                    <Heart className="w-3.5 h-3.5" />
                    {piece.engagement.likes}
                  </div>
                  <div className="flex items-center gap-1 text-xs text-gray-500">
                    <MessageCircle className="w-3.5 h-3.5" />
                    {piece.engagement.comments}
                  </div>
                </div>
              )}
            </div>
          </div>
        ))}
      </div>
    </div>
  );
}
