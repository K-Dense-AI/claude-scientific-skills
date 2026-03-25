"use client";

import { useState } from "react";
import { ChevronLeft, ChevronRight, X } from "lucide-react";

type Post = {
  platform: string;
  type: string;
  caption: string;
  time: string;
};

type ScheduledPosts = Record<string, Post[]>;

const platformColors: Record<string, string> = {
  Instagram: "bg-pink-500",
  Facebook: "bg-blue-500",
  LinkedIn: "bg-cyan-500",
  TikTok: "bg-red-500",
};

const platformDotColors: Record<string, string> = {
  Instagram: "bg-pink-400",
  Facebook: "bg-blue-400",
  LinkedIn: "bg-cyan-400",
  TikTok: "bg-red-400",
};

const dayNames = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"];

function generateMockPosts(): ScheduledPosts {
  const posts: ScheduledPosts = {};
  const now = new Date();
  const year = now.getFullYear();
  const month = now.getMonth();

  const mockEntries: [number, Post[]][] = [
    [2, [
      { platform: "Instagram", type: "Carousel", caption: "Discover our new single-origin beans from Ethiopia. Rich, complex, unforgettable.", time: "09:00" },
      { platform: "Facebook", type: "Post", caption: "Happy Monday! Start your week right with BrewCraft.", time: "10:30" },
    ]],
    [4, [
      { platform: "TikTok", type: "Reel", caption: "POV: You just discovered the perfect pour-over ratio. Game changer.", time: "14:00" },
    ]],
    [5, [
      { platform: "LinkedIn", type: "Post", caption: "How we reduced our carbon footprint by 40% this quarter while scaling production.", time: "08:00" },
    ]],
    [7, [
      { platform: "Instagram", type: "Story", caption: "Behind the scenes at our roasting facility. Swipe up for the full tour!", time: "11:00" },
      { platform: "Instagram", type: "Reel", caption: "3 cold brew recipes you need to try this summer.", time: "16:00" },
    ]],
    [9, [
      { platform: "Facebook", type: "Post", caption: "New blog post: The Science Behind the Perfect Espresso Shot.", time: "09:00" },
      { platform: "TikTok", type: "Reel", caption: "I tried making latte art for 30 days. Here's what happened.", time: "17:00" },
    ]],
    [11, [
      { platform: "Instagram", type: "Post", caption: "Meet Maria, our partner farmer in Colombia. Her family has been growing coffee for 4 generations.", time: "10:00" },
    ]],
    [13, [
      { platform: "Instagram", type: "Carousel", caption: "5 common coffee mistakes and how to fix them. Save this for later!", time: "09:00" },
      { platform: "LinkedIn", type: "Post", caption: "Excited to announce our new partnership with Fair Trade certified cooperatives.", time: "11:00" },
    ]],
    [15, [
      { platform: "TikTok", type: "Reel", caption: "The secret coffee shops don't want you to know about water temperature.", time: "15:00" },
    ]],
    [17, [
      { platform: "Instagram", type: "Story", caption: "Quick poll: Light roast or dark roast? Tell us in the comments!", time: "12:00" },
      { platform: "Facebook", type: "Post", caption: "Weekend vibes. Grab your BrewCraft and relax.", time: "08:00" },
    ]],
    [19, [
      { platform: "Instagram", type: "Reel", caption: "Day in the life of a coffee roaster. It smells even better than you think.", time: "14:00" },
    ]],
    [21, [
      { platform: "Instagram", type: "Post", caption: "Our sustainability report is here. Read about our journey to zero-waste packaging.", time: "09:00" },
      { platform: "LinkedIn", type: "Post", caption: "Q1 community impact: 12 farmer partnerships, 3 new scholarships, 1 mission.", time: "10:00" },
      { platform: "TikTok", type: "Reel", caption: "Stop doing pour-over like this. Do this instead.", time: "16:00" },
    ]],
    [23, [
      { platform: "Instagram", type: "Carousel", caption: "From bean to cup: the complete journey of your morning coffee in 10 slides.", time: "09:30" },
    ]],
    [25, [
      { platform: "Facebook", type: "Post", caption: "Join us this Saturday for our free cupping event at the downtown cafe!", time: "10:00" },
      { platform: "Instagram", type: "Story", caption: "Countdown: 3 days until our new blend drops. Any guesses?", time: "18:00" },
    ]],
    [27, [
      { platform: "Instagram", type: "Post", caption: "The wait is over. Introducing Midnight Harvest - our boldest blend yet.", time: "09:00" },
      { platform: "TikTok", type: "Reel", caption: "First taste reaction of our new Midnight Harvest blend. Was not expecting this.", time: "13:00" },
    ]],
  ];

  mockEntries.forEach(([day, dayPosts]) => {
    const key = `${year}-${String(month + 1).padStart(2, "0")}-${String(day).padStart(2, "0")}`;
    posts[key] = dayPosts;
  });

  return posts;
}

export default function CalendarPage() {
  const [currentDate, setCurrentDate] = useState(new Date());
  const [selectedDay, setSelectedDay] = useState<string | null>(null);
  const scheduledPosts = generateMockPosts();

  const year = currentDate.getFullYear();
  const month = currentDate.getMonth();
  const monthName = currentDate.toLocaleString("default", { month: "long" });

  const firstDay = new Date(year, month, 1);
  const lastDay = new Date(year, month + 1, 0);
  const daysInMonth = lastDay.getDate();

  // Monday = 0, Sunday = 6
  let startDay = firstDay.getDay() - 1;
  if (startDay < 0) startDay = 6;

  const prevMonth = () =>
    setCurrentDate(new Date(year, month - 1, 1));
  const nextMonth = () =>
    setCurrentDate(new Date(year, month + 1, 1));

  const days: (number | null)[] = [];
  for (let i = 0; i < startDay; i++) days.push(null);
  for (let i = 1; i <= daysInMonth; i++) days.push(i);
  while (days.length % 7 !== 0) days.push(null);

  const getPostsForDay = (day: number) => {
    const key = `${year}-${String(month + 1).padStart(2, "0")}-${String(day).padStart(2, "0")}`;
    return scheduledPosts[key] || [];
  };

  const selectedPosts = selectedDay ? scheduledPosts[selectedDay] || [] : [];

  return (
    <div className="p-8">
      <div className="flex items-center justify-between mb-8">
        <div>
          <h1 className="text-3xl font-bold text-white">Content Calendar</h1>
          <p className="text-gray-400 mt-1">
            Plan and schedule your social media content.
          </p>
        </div>
        <div className="flex items-center gap-2">
          {Object.entries(platformDotColors).map(([name, color]) => (
            <div key={name} className="flex items-center gap-1.5 text-xs text-gray-400">
              <div className={`w-2.5 h-2.5 rounded-full ${color}`} />
              {name}
            </div>
          ))}
        </div>
      </div>

      <div className="flex gap-6">
        {/* Calendar Grid */}
        <div className="flex-1 bg-gray-900 border border-gray-800 rounded-2xl p-6">
          {/* Month Navigation */}
          <div className="flex items-center justify-between mb-6">
            <button
              onClick={prevMonth}
              className="p-2 hover:bg-gray-800 rounded-xl transition-colors"
            >
              <ChevronLeft className="w-5 h-5 text-gray-400" />
            </button>
            <h2 className="text-lg font-semibold text-white">
              {monthName} {year}
            </h2>
            <button
              onClick={nextMonth}
              className="p-2 hover:bg-gray-800 rounded-xl transition-colors"
            >
              <ChevronRight className="w-5 h-5 text-gray-400" />
            </button>
          </div>

          {/* Day Headers */}
          <div className="grid grid-cols-7 gap-1 mb-2">
            {dayNames.map((d) => (
              <div
                key={d}
                className="text-center text-xs font-medium text-gray-500 py-2"
              >
                {d}
              </div>
            ))}
          </div>

          {/* Day Cells */}
          <div className="grid grid-cols-7 gap-1">
            {days.map((day, i) => {
              if (day === null) {
                return <div key={`empty-${i}`} className="h-24" />;
              }
              const posts = getPostsForDay(day);
              const dayKey = `${year}-${String(month + 1).padStart(2, "0")}-${String(day).padStart(2, "0")}`;
              const isSelected = selectedDay === dayKey;
              const isToday =
                day === new Date().getDate() &&
                month === new Date().getMonth() &&
                year === new Date().getFullYear();

              return (
                <button
                  key={`day-${day}`}
                  onClick={() => setSelectedDay(isSelected ? null : dayKey)}
                  className={`h-24 p-2 rounded-xl text-left transition-all border ${
                    isSelected
                      ? "border-purple-500 bg-purple-500/5"
                      : "border-transparent hover:bg-gray-800/50"
                  }`}
                >
                  <span
                    className={`text-sm font-medium ${
                      isToday
                        ? "bg-purple-600 text-white w-6 h-6 rounded-full flex items-center justify-center"
                        : "text-gray-300"
                    }`}
                  >
                    {day}
                  </span>
                  {posts.length > 0 && (
                    <div className="flex flex-wrap gap-1 mt-2">
                      {posts.map((post, j) => (
                        <div
                          key={j}
                          className={`w-2 h-2 rounded-full ${platformDotColors[post.platform]}`}
                        />
                      ))}
                    </div>
                  )}
                </button>
              );
            })}
          </div>
        </div>

        {/* Side Panel */}
        {selectedDay && (
          <div className="w-80 bg-gray-900 border border-gray-800 rounded-2xl p-5 h-fit">
            <div className="flex items-center justify-between mb-4">
              <h3 className="text-sm font-semibold text-white">
                {new Date(selectedDay + "T00:00:00").toLocaleDateString(
                  "default",
                  {
                    weekday: "long",
                    month: "long",
                    day: "numeric",
                  }
                )}
              </h3>
              <button
                onClick={() => setSelectedDay(null)}
                className="text-gray-500 hover:text-white"
              >
                <X className="w-4 h-4" />
              </button>
            </div>

            {selectedPosts.length === 0 ? (
              <p className="text-sm text-gray-500">
                No posts scheduled for this day.
              </p>
            ) : (
              <div className="space-y-3">
                {selectedPosts.map((post, i) => (
                  <div
                    key={i}
                    className="bg-gray-800/50 rounded-xl p-4 space-y-2"
                  >
                    <div className="flex items-center gap-2">
                      <span
                        className={`text-xs px-2 py-0.5 rounded-full text-white font-medium ${platformColors[post.platform]}`}
                      >
                        {post.platform}
                      </span>
                      <span className="text-xs px-2 py-0.5 rounded-full bg-gray-700 text-gray-300">
                        {post.type}
                      </span>
                    </div>
                    <p className="text-sm text-gray-300 leading-relaxed">
                      {post.caption}
                    </p>
                    <p className="text-xs text-gray-500">{post.time}</p>
                  </div>
                ))}
              </div>
            )}
          </div>
        )}
      </div>
    </div>
  );
}
