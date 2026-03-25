import { NextRequest, NextResponse } from "next/server";

export async function POST(request: NextRequest) {
  const body = await request.json();

  // Simulate processing delay
  await new Promise((resolve) => setTimeout(resolve, 2000));

  const campaignPlan = {
    id: `campaign_${Date.now()}`,
    name: body.name || "Untitled Campaign",
    type: body.type || "brand_awareness",
    platforms: body.platforms || ["instagram"],
    duration: body.duration || 30,
    status: "generated",
    contentPlan: {
      totalPosts: Math.floor((body.duration || 30) * 0.8),
      breakdown: {
        posts: Math.floor((body.duration || 30) * 0.3),
        reels: Math.floor((body.duration || 30) * 0.25),
        stories: Math.floor((body.duration || 30) * 0.2),
        carousels: Math.floor((body.duration || 30) * 0.05),
      },
    },
    contentPillars: body.pillars || ["Behind the Scenes", "Education", "Community"],
    hooks: body.hooks || [],
    tone: body.tone || "professional",
    createdAt: new Date().toISOString(),
  };

  return NextResponse.json(campaignPlan);
}
