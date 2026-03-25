import { NextRequest, NextResponse } from "next/server";

interface BrandDna {
  url: string;
  name: string;
  tagline: string;
  industry: string;
  colors: { primary: string; secondary: string; accent: string; neutral: string };
  typography: { heading: string; body: string };
  voice: { traits: string[]; tone: string; doNots: string[] };
  audience: { primary: string; interests: string[] };
  contentPillars: string[];
  socialProfiles: Record<string, string>;
  analyzedAt: string;
}

function extractColors(html: string): string[] {
  const colorRegex = /#(?:[0-9a-fA-F]{6}|[0-9a-fA-F]{3})\b/g;
  const matches = html.match(colorRegex) || [];
  const unique = [...new Set(matches.map((c) => c.toUpperCase()))];
  // Filter out common white/black/gray
  const filtered = unique.filter(
    (c) => !["#FFFFFF", "#FFF", "#000000", "#000", "#333333", "#333", "#666666", "#666", "#999999", "#999", "#CCCCCC", "#CCC"].includes(c)
  );
  return filtered.slice(0, 6);
}

function extractTitle(html: string, url: string): string {
  // Prefer og:site_name (most reliable brand name)
  const ogSite = html.match(/property="og:site_name"\s+content="([^"]+)"/i)
    || html.match(/content="([^"]+)"\s+property="og:site_name"/i);
  if (ogSite) return ogSite[1].trim();

  // Try og:title
  const ogTitle = html.match(/property="og:title"\s+content="([^"]+)"/i)
    || html.match(/content="([^"]+)"\s+property="og:title"/i);
  if (ogTitle) {
    const t = ogTitle[1].trim();
    if (!["acasa", "home", "homepage", "index"].includes(t.toLowerCase())) return t;
  }

  // Try <title> but skip generic titles
  const titleMatch = html.match(/<title[^>]*>([^<]+)<\/title>/i);
  if (titleMatch) {
    const t = titleMatch[1].trim().split(/[|\-–—]/)[0].trim();
    if (!["acasa", "home", "homepage", "index"].includes(t.toLowerCase())) return t;
  }

  // Fallback: extract brand from URL
  try {
    const hostname = new URL(url.startsWith("http") ? url : "https://" + url).hostname.replace("www.", "");
    const name = hostname.split(".")[0];
    return name.charAt(0).toUpperCase() + name.slice(1);
  } catch {
    return "";
  }
}

function extractDescription(html: string): string {
  const metaDesc = html.match(/name="description"\s+content="([^"]+)"/i)
    || html.match(/property="og:description"\s+content="([^"]+)"/i);
  return metaDesc ? metaDesc[1].trim() : "";
}

function extractFonts(html: string): { heading: string; body: string } {
  const googleFonts = html.match(/fonts\.googleapis\.com\/css2?\?family=([^"&]+)/gi) || [];
  const fontNames: string[] = [];
  for (const url of googleFonts) {
    const fam = url.match(/family=([^:&"]+)/);
    if (fam) fontNames.push(decodeURIComponent(fam[1]).replace(/\+/g, " "));
  }
  const cssFonts = html.match(/font-family:\s*['"]?([^'";,}]+)/gi) || [];
  for (const f of cssFonts) {
    const name = f.replace(/font-family:\s*['"]?/i, "").trim();
    if (name && !name.includes("inherit") && !name.includes("sans-serif") && !name.includes("serif") && name.length < 30) {
      fontNames.push(name);
    }
  }
  const unique = [...new Set(fontNames)];
  return {
    heading: unique[0] || "System Default",
    body: unique[1] || unique[0] || "System Default",
  };
}

function guessIndustry(text: string): string {
  const lower = text.toLowerCase();
  const industries: [string, string[]][] = [
    ["Pet Care & Pet Food", ["pet", "dog", "cat", "puppy", "kitten", "animal", "catel", "caine", "pisica", "hrana", "mancare de", "food for dogs", "food for cats", "pet food", "pet shop"]],
    ["Food & Beverage", ["food", "restaurant", "coffee", "tea", "cafe", "drink", "cuisine", "recipe", "cook", "bakery", "pizza", "burger"]],
    ["Technology", ["software", "app", "tech", "digital", "cloud", "saas", "platform", "ai", "machine learning"]],
    ["Health & Wellness", ["health", "wellness", "fitness", "gym", "yoga", "meditation", "supplement", "vitamin"]],
    ["Fashion & Beauty", ["fashion", "beauty", "clothing", "apparel", "cosmetic", "skincare", "makeup", "style"]],
    ["E-commerce & Retail", ["shop", "store", "buy", "cart", "product", "order", "delivery", "shipping"]],
    ["Education", ["learn", "course", "education", "training", "academy", "school", "university", "student"]],
    ["Finance", ["finance", "bank", "invest", "insurance", "loan", "credit", "money", "fintech"]],
    ["Real Estate", ["property", "real estate", "apartment", "house", "rent", "home", "imobiliare"]],
    ["Travel & Tourism", ["travel", "tourism", "hotel", "flight", "vacation", "booking", "trip"]],
    ["Marketing & Advertising", ["marketing", "advertising", "brand", "agency", "creative", "campaign"]],
  ];
  let bestMatch = "General Business";
  let bestScore = 0;
  for (const [industry, keywords] of industries) {
    let score = 0;
    for (const kw of keywords) {
      if (lower.includes(kw)) score++;
    }
    if (score > bestScore) {
      bestScore = score;
      bestMatch = industry;
    }
  }
  return bestMatch;
}

function guessVoiceTraits(text: string): string[] {
  const lower = text.toLowerCase();
  const traits: string[] = [];
  if (lower.includes("premium") || lower.includes("luxury") || lower.includes("exclusive")) traits.push("Premium");
  if (lower.includes("natural") || lower.includes("organic") || lower.includes("eco")) traits.push("Eco-Conscious");
  if (lower.includes("fun") || lower.includes("play") || lower.includes("joy") || lower.includes("happy")) traits.push("Playful");
  if (lower.includes("expert") || lower.includes("professional") || lower.includes("specialist")) traits.push("Expert");
  if (lower.includes("love") || lower.includes("care") || lower.includes("passion") || lower.includes("heart") || lower.includes("iubire") || lower.includes("grija")) traits.push("Caring");
  if (lower.includes("innovat") || lower.includes("modern") || lower.includes("new")) traits.push("Innovative");
  if (lower.includes("trust") || lower.includes("reliable") || lower.includes("quality") || lower.includes("calitate")) traits.push("Trustworthy");
  if (lower.includes("community") || lower.includes("together") || lower.includes("join") || lower.includes("comunitate")) traits.push("Community-Oriented");
  if (traits.length === 0) traits.push("Warm", "Authentic");
  return traits.slice(0, 4);
}

function extractContentPillars(text: string, industry: string): string[] {
  const lower = text.toLowerCase();
  const pillars: string[] = [];

  if (industry.includes("Pet")) {
    if (lower.includes("nutrition") || lower.includes("hrana") || lower.includes("food") || lower.includes("mancare")) pillars.push("Nutritie & Sanatate");
    if (lower.includes("recipe") || lower.includes("reteta") || lower.includes("ingredient")) pillars.push("Retete & Ingrediente");
    if (lower.includes("tips") || lower.includes("sfat") || lower.includes("advice") || lower.includes("care")) pillars.push("Sfaturi & Ingrijire");
    if (lower.includes("story") || lower.includes("povest") || lower.includes("happy") || lower.includes("love")) pillars.push("Povesti de Succes");
    if (pillars.length < 3) pillars.push("Behind the Brand");
  } else {
    if (lower.includes("about") || lower.includes("story") || lower.includes("mission")) pillars.push("Brand Story");
    if (lower.includes("product") || lower.includes("service") || lower.includes("feature")) pillars.push("Products & Services");
    if (lower.includes("blog") || lower.includes("news") || lower.includes("article")) pillars.push("Educational Content");
    if (lower.includes("team") || lower.includes("people") || lower.includes("culture")) pillars.push("Behind the Scenes");
    if (lower.includes("testimonial") || lower.includes("review") || lower.includes("client")) pillars.push("Social Proof");
  }

  if (pillars.length === 0) pillars.push("Brand Story", "Products & Services", "Educational Content");
  return pillars.slice(0, 5);
}

export async function POST(request: NextRequest) {
  const body = await request.json();
  const { url, socialHandles } = body;

  if (!url) {
    return NextResponse.json({ error: "URL is required" }, { status: 400 });
  }

  // Normalize URL
  let normalizedUrl = url.trim();
  if (!normalizedUrl.startsWith("http")) {
    normalizedUrl = "https://" + normalizedUrl;
  }

  try {
    // Fetch the actual website
    const response = await fetch(normalizedUrl, {
      headers: {
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        "Accept-Language": "en-US,en;q=0.9,ro;q=0.8",
      },
      signal: AbortSignal.timeout(20000),
      redirect: "follow",
    });

    if (!response.ok && response.status !== 503) {
      throw new Error(`Failed to fetch: ${response.status}`);
    }
    // Some sites return 503 with content (Cloudflare challenge pages with actual HTML behind)
    // Try to read content anyway

    const html = await response.text();

    // Extract brand information from HTML
    const rawName = extractTitle(html, normalizedUrl);
    const description = extractDescription(html);
    const colors = extractColors(html);
    const fonts = extractFonts(html);
    const fullText = html.replace(/<[^>]+>/g, " ").replace(/\s+/g, " ");
    const industry = guessIndustry(fullText + " " + description);
    const voiceTraits = guessVoiceTraits(fullText + " " + description);
    const contentPillars = extractContentPillars(fullText, industry);

    // Build color palette
    const colorPalette = {
      primary: colors[0] || "#6B21A8",
      secondary: colors[1] || "#3B82F6",
      accent: colors[2] || "#10B981",
      neutral: colors[3] || "#F3F4F6",
    };

    const brandDna: BrandDna = {
      url: normalizedUrl,
      name: rawName || new URL(normalizedUrl).hostname.replace("www.", "").split(".")[0],
      tagline: description.slice(0, 120) || "No tagline detected",
      industry,
      colors: colorPalette,
      typography: fonts,
      voice: {
        traits: voiceTraits,
        tone: voiceTraits.includes("Playful") ? "Casual and fun" :
              voiceTraits.includes("Expert") ? "Authoritative yet approachable" :
              voiceTraits.includes("Premium") ? "Elegant and refined" :
              "Warm and authentic",
        doNots: [
          "Avoid generic corporate language",
          "Don't use competitors' slogans",
          "No misleading claims",
        ],
      },
      audience: {
        primary: industry.includes("Pet") ? "Pet owners, 25-55, urban & suburban" :
                 industry.includes("Tech") ? "Tech professionals, 25-45, early adopters" :
                 industry.includes("Fashion") ? "Style-conscious consumers, 18-40" :
                 "General consumers, 25-55",
        interests: contentPillars.slice(0, 4),
      },
      contentPillars,
      socialProfiles: socialHandles || {},
      analyzedAt: new Date().toISOString(),
    };

    return NextResponse.json(brandDna);
  } catch (error) {
    // If fetch fails, try to build a basic profile from the URL alone
    const hostname = (() => {
      try { return new URL(normalizedUrl).hostname.replace("www.", ""); }
      catch { return normalizedUrl; }
    })();
    const brandName = hostname.split(".")[0];
    const capitalized = brandName.charAt(0).toUpperCase() + brandName.slice(1);

    const fallbackDna: BrandDna = {
      url: normalizedUrl,
      name: capitalized,
      tagline: `Welcome to ${capitalized}`,
      industry: "General Business",
      colors: {
        primary: "#6B21A8",
        secondary: "#3B82F6",
        accent: "#10B981",
        neutral: "#F3F4F6",
      },
      typography: { heading: "System Default", body: "System Default" },
      voice: {
        traits: ["Professional", "Approachable"],
        tone: "Professional yet friendly",
        doNots: ["Avoid jargon", "No aggressive sales language"],
      },
      audience: {
        primary: "General consumers, 25-55",
        interests: ["Brand Story", "Products", "Education"],
      },
      contentPillars: ["Brand Story", "Products & Services", "Educational Content"],
      socialProfiles: socialHandles || {},
      analyzedAt: new Date().toISOString(),
    };

    return NextResponse.json(fallbackDna);
  }
}
