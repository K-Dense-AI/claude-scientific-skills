"use client";

import { useState } from "react";
import { usePathname } from "next/navigation";
import Link from "next/link";
import {
  Sparkles,
  LayoutDashboard,
  PlusCircle,
  Dna,
  Calendar,
  Library,
  BarChart3,
  Megaphone,
  Settings,
  ChevronLeft,
  ChevronRight,
} from "lucide-react";

const navItems = [
  { label: "Dashboard", href: "/", icon: LayoutDashboard },
  { label: "New Campaign", href: "/campaign/new", icon: PlusCircle },
  { label: "Brand DNA", href: "/brand", icon: Dna },
  { label: "Content Calendar", href: "/calendar", icon: Calendar },
  { label: "Content Library", href: "/content", icon: Library },
  { label: "Ads & Publishing", href: "/ads", icon: Megaphone },
  { label: "Analytics", href: "/analytics", icon: BarChart3 },
  { label: "Settings", href: "/settings", icon: Settings },
];

export default function Sidebar() {
  const [collapsed, setCollapsed] = useState(false);
  const pathname = usePathname();

  return (
    <aside
      className={`flex flex-col bg-gray-900 border-r border-gray-800 transition-all duration-300 ${
        collapsed ? "w-20" : "w-64"
      } min-h-screen`}
    >
      {/* Logo */}
      <div className="flex items-center gap-3 px-5 py-6 border-b border-gray-800">
        <div className="flex items-center justify-center w-10 h-10 rounded-xl bg-purple-600">
          <Sparkles className="w-5 h-5 text-white" />
        </div>
        {!collapsed && (
          <span className="text-xl font-bold text-white tracking-tight">
            SocialAI
          </span>
        )}
      </div>

      {/* Navigation */}
      <nav className="flex-1 px-3 py-4 space-y-1">
        {navItems.map((item) => {
          const isActive =
            pathname === item.href ||
            (item.href !== "/" && pathname.startsWith(item.href));
          const Icon = item.icon;

          return (
            <Link
              key={item.href}
              href={item.href}
              className={`flex items-center gap-3 px-3 py-2.5 rounded-xl text-sm font-medium transition-all duration-200 ${
                isActive
                  ? "bg-purple-600/20 text-purple-400 border border-purple-500/30"
                  : "text-gray-400 hover:text-white hover:bg-gray-800"
              } ${collapsed ? "justify-center" : ""}`}
            >
              <Icon className="w-5 h-5 shrink-0" />
              {!collapsed && <span>{item.label}</span>}
            </Link>
          );
        })}
      </nav>

      {/* Collapse toggle */}
      <button
        onClick={() => setCollapsed(!collapsed)}
        className="flex items-center justify-center gap-2 px-3 py-4 text-gray-400 hover:text-white border-t border-gray-800 transition-colors"
      >
        {collapsed ? (
          <ChevronRight className="w-5 h-5" />
        ) : (
          <>
            <ChevronLeft className="w-5 h-5" />
            <span className="text-sm">Collapse</span>
          </>
        )}
      </button>
    </aside>
  );
}
