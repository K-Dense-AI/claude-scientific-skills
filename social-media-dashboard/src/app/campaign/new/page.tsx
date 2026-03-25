"use client";

import { useState } from "react";
import {
  Globe,
  Instagram,
  Facebook,
  Linkedin,
  ArrowRight,
  ArrowLeft,
  Check,
  Loader2,
  Sparkles,
  Calendar,
  Library,
} from "lucide-react";
import Link from "next/link";
import {
  PieChart,
  Pie,
  Cell,
  ResponsiveContainer,
  Tooltip,
} from "recharts";

const steps = [
  "Brand Discovery",
  "Campaign Setup",
  "Content Strategy",
  "Review & Generate",
];

const campaignTypes = [
  "Product Launch",
  "Brand Awareness",
  "Engagement Growth",
  "Lead Generation",
  "Event Promotion",
  "Seasonal",
];

const platformOptions = [
  { id: "instagram", label: "Instagram", color: "bg-pink-500" },
  { id: "facebook", label: "Facebook", color: "bg-blue-500" },
  { id: "linkedin", label: "LinkedIn", color: "bg-cyan-500" },
  { id: "tiktok", label: "TikTok", color: "bg-red-500" },
];

const viralHooks = [
  "How to [achieve X] in [time]",
  "I tried [X] for 30 days",
  "The secret [industry] doesn't want you to know",
  "Stop doing [X] — do this instead",
  "POV: You just discovered [X]",
  "[Number] things I wish I knew about [X]",
  "Why [popular opinion] is wrong",
  "The [X] hack that changed my life",
  "Day in the life of [persona]",
  "What [X] looks like vs what it actually is",
];

const toneOptions = ["Professional", "Casual", "Bold", "Playful", "Inspirational"];

const contentMixData = [
  { name: "Value Content", value: 60, color: "#a855f7" },
  { name: "Engagement Content", value: 30, color: "#3b82f6" },
  { name: "Promotional Content", value: 10, color: "#f43f5e" },
];

const generationSteps = [
  "Mining brand identity...",
  "Creating content briefs...",
  "Generating visuals...",
  "Building calendar...",
];

export default function NewCampaignPage() {
  const [step, setStep] = useState(0);

  // Step 1 state
  const [brandUrl, setBrandUrl] = useState("");
  const [socialHandles, setSocialHandles] = useState({
    instagram: "",
    facebook: "",
    linkedin: "",
    tiktok: "",
  });
  const [analyzingBrand, setAnalyzingBrand] = useState(false);
  const [brandAnalyzed, setBrandAnalyzed] = useState(false);
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const [brandData, setBrandData] = useState<any>(null);

  // Step 2 state
  const [campaignName, setCampaignName] = useState("");
  const [campaignType, setCampaignType] = useState("Brand Awareness");
  const [selectedPlatforms, setSelectedPlatforms] = useState<string[]>([
    "instagram",
  ]);
  const [duration, setDuration] = useState(30);
  const [budget, setBudget] = useState("");

  // Step 3 state
  const [pillars, setPillars] = useState<string[]>([
    "Behind the Scenes",
    "Education",
  ]);
  const [pillarInput, setPillarInput] = useState("");
  const [selectedHooks, setSelectedHooks] = useState<string[]>([]);
  const [tone, setTone] = useState("Professional");

  // Step 4 state
  const [generating, setGenerating] = useState(false);
  const [genStep, setGenStep] = useState(0);
  const [generated, setGenerated] = useState(false);

  const analyzeBrand = async () => {
    setAnalyzingBrand(true);
    try {
      const res = await fetch("/api/brand", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ url: brandUrl, socialHandles }),
      });
      const data = await res.json();
      setBrandData(data);
      setBrandAnalyzed(true);
      // Auto-fill pillars from brand analysis
      if (data.contentPillars?.length) {
        setPillars(data.contentPillars);
      }
    } catch {
      setBrandAnalyzed(true);
      setBrandData({ name: "Unknown", industry: "General", voice: { traits: ["Professional"] }, colors: { primary: "#6B21A8", secondary: "#3B82F6", accent: "#10B981", neutral: "#F3F4F6" } });
    }
    setAnalyzingBrand(false);
  };

  const addPillar = () => {
    if (pillarInput.trim() && !pillars.includes(pillarInput.trim())) {
      setPillars([...pillars, pillarInput.trim()]);
      setPillarInput("");
    }
  };

  const removePillar = (p: string) => {
    setPillars(pillars.filter((x) => x !== p));
  };

  const togglePlatform = (id: string) => {
    setSelectedPlatforms((prev) =>
      prev.includes(id) ? prev.filter((p) => p !== id) : [...prev, id]
    );
  };

  const toggleHook = (hook: string) => {
    setSelectedHooks((prev) =>
      prev.includes(hook) ? prev.filter((h) => h !== hook) : [...prev, hook]
    );
  };

  const generateCampaign = async () => {
    setGenerating(true);
    for (let i = 0; i < generationSteps.length; i++) {
      setGenStep(i);
      await new Promise((r) => setTimeout(r, 1500));
    }
    setGenerating(false);
    setGenerated(true);
  };

  return (
    <div className="p-8 max-w-4xl mx-auto">
      <h1 className="text-3xl font-bold text-white mb-2">New Campaign</h1>
      <p className="text-gray-400 mb-8">
        Create an AI-powered social media campaign in minutes.
      </p>

      {/* Step Indicator */}
      <div className="flex items-center gap-2 mb-10">
        {steps.map((s, i) => (
          <div key={s} className="flex items-center gap-2 flex-1">
            <div
              className={`flex items-center justify-center w-8 h-8 rounded-full text-sm font-bold shrink-0 transition-colors ${
                i < step
                  ? "bg-purple-600 text-white"
                  : i === step
                  ? "bg-purple-500 text-white"
                  : "bg-gray-800 text-gray-500"
              }`}
            >
              {i < step ? <Check className="w-4 h-4" /> : i + 1}
            </div>
            <span
              className={`text-sm hidden sm:block ${
                i === step ? "text-white font-medium" : "text-gray-500"
              }`}
            >
              {s}
            </span>
            {i < steps.length - 1 && (
              <div
                className={`flex-1 h-px ${
                  i < step ? "bg-purple-600" : "bg-gray-800"
                }`}
              />
            )}
          </div>
        ))}
      </div>

      {/* Step 1: Brand Discovery */}
      {step === 0 && (
        <div className="space-y-6 transition-all duration-300">
          <div className="bg-gray-900 border border-gray-800 rounded-2xl p-6 space-y-5">
            <h2 className="text-xl font-semibold text-white">
              Brand Discovery
            </h2>
            <p className="text-sm text-gray-400">
              Enter your brand&apos;s website and social handles so our AI can
              analyze your brand identity.
            </p>

            <div>
              <label className="block text-sm font-medium text-gray-300 mb-1.5">
                Brand Website URL
              </label>
              <div className="flex items-center gap-2 bg-gray-800 rounded-xl px-4 py-3 border border-gray-700 focus-within:border-purple-500">
                <Globe className="w-4 h-4 text-gray-500" />
                <input
                  type="url"
                  value={brandUrl}
                  onChange={(e) => setBrandUrl(e.target.value)}
                  placeholder="https://yourwebsite.com"
                  className="bg-transparent flex-1 text-white text-sm outline-none placeholder:text-gray-600"
                />
              </div>
            </div>

            <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
              {[
                { key: "instagram", label: "Instagram", Icon: Instagram },
                { key: "facebook", label: "Facebook", Icon: Facebook },
                { key: "linkedin", label: "LinkedIn", Icon: Linkedin },
                { key: "tiktok", label: "TikTok", Icon: Globe },
              ].map(({ key, label, Icon }) => (
                <div key={key}>
                  <label className="block text-sm font-medium text-gray-300 mb-1.5">
                    {label}
                  </label>
                  <div className="flex items-center gap-2 bg-gray-800 rounded-xl px-4 py-3 border border-gray-700 focus-within:border-purple-500">
                    <Icon className="w-4 h-4 text-gray-500" />
                    <input
                      type="text"
                      value={socialHandles[key as keyof typeof socialHandles]}
                      onChange={(e) =>
                        setSocialHandles({
                          ...socialHandles,
                          [key]: e.target.value,
                        })
                      }
                      placeholder={`@${label.toLowerCase()}_handle`}
                      className="bg-transparent flex-1 text-white text-sm outline-none placeholder:text-gray-600"
                    />
                  </div>
                </div>
              ))}
            </div>

            <button
              onClick={analyzeBrand}
              disabled={analyzingBrand || !brandUrl}
              className="w-full py-3 px-6 bg-purple-600 hover:bg-purple-700 disabled:bg-gray-700 disabled:text-gray-500 text-white font-medium rounded-xl transition-colors flex items-center justify-center gap-2"
            >
              {analyzingBrand ? (
                <>
                  <Loader2 className="w-4 h-4 animate-spin" />
                  Analyzing Brand...
                </>
              ) : (
                <>
                  <Sparkles className="w-4 h-4" />
                  Analyze Brand
                </>
              )}
            </button>

            {brandAnalyzed && brandData && (
              <div className="bg-gray-800/50 rounded-xl p-5 border border-green-500/20 space-y-3">
                <div className="flex items-center gap-2 text-green-400 text-sm font-medium">
                  <Check className="w-4 h-4" />
                  Brand DNA extracted successfully
                </div>
                <div className="grid grid-cols-2 gap-4 text-sm">
                  <div>
                    <p className="text-gray-500">Brand Name</p>
                    <p className="text-white font-medium">{brandData.name}</p>
                  </div>
                  <div>
                    <p className="text-gray-500">Industry</p>
                    <p className="text-white">{brandData.industry}</p>
                  </div>
                  <div>
                    <p className="text-gray-500">Voice</p>
                    <p className="text-white">{brandData.voice?.traits?.join(", ") || "N/A"}</p>
                  </div>
                  <div>
                    <p className="text-gray-500">Colors Found</p>
                    <div className="flex gap-1.5 mt-1">
                      {Object.values(brandData.colors || {}).map(
                        (c) => (
                          <div
                            key={c as string}
                            className="w-5 h-5 rounded-full border border-gray-600"
                            style={{ backgroundColor: c as string }}
                            title={c as string}
                          />
                        )
                      )}
                    </div>
                  </div>
                  {brandData.tagline && brandData.tagline !== "No tagline detected" && (
                    <div className="col-span-2">
                      <p className="text-gray-500">Tagline</p>
                      <p className="text-white text-xs">{brandData.tagline}</p>
                    </div>
                  )}
                  {brandData.contentPillars?.length > 0 && (
                    <div className="col-span-2">
                      <p className="text-gray-500 mb-1">Content Pillars (auto-detected)</p>
                      <div className="flex flex-wrap gap-1.5">
                        {brandData.contentPillars.map((p: string) => (
                          <span key={p} className="px-2 py-0.5 bg-purple-500/10 text-purple-400 rounded-full text-xs">{p}</span>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
              </div>
            )}
          </div>
        </div>
      )}

      {/* Step 2: Campaign Setup */}
      {step === 1 && (
        <div className="space-y-6 transition-all duration-300">
          <div className="bg-gray-900 border border-gray-800 rounded-2xl p-6 space-y-5">
            <h2 className="text-xl font-semibold text-white">
              Campaign Setup
            </h2>

            <div>
              <label className="block text-sm font-medium text-gray-300 mb-1.5">
                Campaign Name
              </label>
              <input
                type="text"
                value={campaignName}
                onChange={(e) => setCampaignName(e.target.value)}
                placeholder="e.g., Summer Launch 2025"
                className="w-full bg-gray-800 rounded-xl px-4 py-3 border border-gray-700 text-white text-sm outline-none focus:border-purple-500 placeholder:text-gray-600"
              />
            </div>

            <div>
              <label className="block text-sm font-medium text-gray-300 mb-1.5">
                Campaign Type
              </label>
              <select
                value={campaignType}
                onChange={(e) => setCampaignType(e.target.value)}
                className="w-full bg-gray-800 rounded-xl px-4 py-3 border border-gray-700 text-white text-sm outline-none focus:border-purple-500"
              >
                {campaignTypes.map((t) => (
                  <option key={t} value={t}>
                    {t}
                  </option>
                ))}
              </select>
            </div>

            <div>
              <label className="block text-sm font-medium text-gray-300 mb-2">
                Target Platforms
              </label>
              <div className="flex flex-wrap gap-3">
                {platformOptions.map((p) => (
                  <button
                    key={p.id}
                    onClick={() => togglePlatform(p.id)}
                    className={`flex items-center gap-2 px-4 py-2.5 rounded-xl text-sm font-medium transition-all border ${
                      selectedPlatforms.includes(p.id)
                        ? "border-purple-500 bg-purple-500/10 text-purple-300"
                        : "border-gray-700 bg-gray-800 text-gray-400 hover:border-gray-600"
                    }`}
                  >
                    <div className={`w-3 h-3 rounded-full ${p.color}`} />
                    {p.label}
                  </button>
                ))}
              </div>
            </div>

            <div>
              <label className="block text-sm font-medium text-gray-300 mb-1.5">
                Duration: {duration} days
              </label>
              <input
                type="range"
                min={7}
                max={90}
                value={duration}
                onChange={(e) => setDuration(Number(e.target.value))}
                className="w-full accent-purple-500"
              />
              <div className="flex justify-between text-xs text-gray-500 mt-1">
                <span>7 days</span>
                <span>90 days</span>
              </div>
            </div>

            <div>
              <label className="block text-sm font-medium text-gray-300 mb-1.5">
                Budget (optional)
              </label>
              <input
                type="text"
                value={budget}
                onChange={(e) => setBudget(e.target.value)}
                placeholder="$5,000"
                className="w-full bg-gray-800 rounded-xl px-4 py-3 border border-gray-700 text-white text-sm outline-none focus:border-purple-500 placeholder:text-gray-600"
              />
            </div>
          </div>
        </div>
      )}

      {/* Step 3: Content Strategy */}
      {step === 2 && (
        <div className="space-y-6 transition-all duration-300">
          <div className="bg-gray-900 border border-gray-800 rounded-2xl p-6 space-y-5">
            <h2 className="text-xl font-semibold text-white">
              Content Strategy
            </h2>

            {/* Content Pillars */}
            <div>
              <label className="block text-sm font-medium text-gray-300 mb-2">
                Content Pillars
              </label>
              <div className="flex flex-wrap gap-2 mb-3">
                {pillars.map((p) => (
                  <span
                    key={p}
                    className="flex items-center gap-1.5 px-3 py-1.5 bg-purple-500/10 text-purple-400 rounded-full text-sm"
                  >
                    {p}
                    <button
                      onClick={() => removePillar(p)}
                      className="hover:text-purple-200"
                    >
                      x
                    </button>
                  </span>
                ))}
              </div>
              <div className="flex gap-2">
                <input
                  type="text"
                  value={pillarInput}
                  onChange={(e) => setPillarInput(e.target.value)}
                  onKeyDown={(e) => e.key === "Enter" && addPillar()}
                  placeholder="Add a content pillar..."
                  className="flex-1 bg-gray-800 rounded-xl px-4 py-2.5 border border-gray-700 text-white text-sm outline-none focus:border-purple-500 placeholder:text-gray-600"
                />
                <button
                  onClick={addPillar}
                  className="px-4 py-2.5 bg-gray-800 hover:bg-gray-700 text-white rounded-xl text-sm border border-gray-700"
                >
                  Add
                </button>
              </div>
            </div>

            {/* Content Mix */}
            <div>
              <label className="block text-sm font-medium text-gray-300 mb-2">
                Content Mix (60/30/10 Rule)
              </label>
              <div className="flex items-center gap-6">
                <div className="w-40 h-40">
                  <ResponsiveContainer width="100%" height="100%">
                    <PieChart>
                      <Pie
                        data={contentMixData}
                        cx="50%"
                        cy="50%"
                        innerRadius={35}
                        outerRadius={65}
                        dataKey="value"
                        stroke="none"
                      >
                        {contentMixData.map((entry, index) => (
                          <Cell key={index} fill={entry.color} />
                        ))}
                      </Pie>
                      <Tooltip
                        contentStyle={{
                          backgroundColor: "#1f2937",
                          border: "1px solid #374151",
                          borderRadius: "8px",
                          color: "#fff",
                          fontSize: "12px",
                        }}
                      />
                    </PieChart>
                  </ResponsiveContainer>
                </div>
                <div className="space-y-2">
                  {contentMixData.map((item) => (
                    <div key={item.name} className="flex items-center gap-2">
                      <div
                        className="w-3 h-3 rounded-full"
                        style={{ backgroundColor: item.color }}
                      />
                      <span className="text-sm text-gray-300">
                        {item.name} — {item.value}%
                      </span>
                    </div>
                  ))}
                </div>
              </div>
            </div>

            {/* Viral Hooks */}
            <div>
              <label className="block text-sm font-medium text-gray-300 mb-2">
                Viral Hooks
              </label>
              <div className="grid grid-cols-1 sm:grid-cols-2 gap-2">
                {viralHooks.map((hook) => (
                  <button
                    key={hook}
                    onClick={() => toggleHook(hook)}
                    className={`text-left px-3 py-2.5 rounded-xl text-sm transition-all border ${
                      selectedHooks.includes(hook)
                        ? "border-purple-500 bg-purple-500/10 text-purple-300"
                        : "border-gray-700 bg-gray-800 text-gray-400 hover:border-gray-600"
                    }`}
                  >
                    {hook}
                  </button>
                ))}
              </div>
            </div>

            {/* Tone */}
            <div>
              <label className="block text-sm font-medium text-gray-300 mb-2">
                Tone of Voice
              </label>
              <div className="flex flex-wrap gap-2">
                {toneOptions.map((t) => (
                  <button
                    key={t}
                    onClick={() => setTone(t)}
                    className={`px-4 py-2 rounded-xl text-sm font-medium transition-all border ${
                      tone === t
                        ? "border-purple-500 bg-purple-600 text-white"
                        : "border-gray-700 bg-gray-800 text-gray-400 hover:border-gray-600"
                    }`}
                  >
                    {t}
                  </button>
                ))}
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Step 4: Review & Generate */}
      {step === 3 && (
        <div className="space-y-6 transition-all duration-300">
          {!generated ? (
            <div className="bg-gray-900 border border-gray-800 rounded-2xl p-6 space-y-5">
              <h2 className="text-xl font-semibold text-white">
                Review &amp; Generate
              </h2>

              {/* Summary */}
              <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
                <div className="bg-gray-800/50 rounded-xl p-4">
                  <p className="text-xs text-gray-500 uppercase tracking-wider mb-1">
                    Campaign
                  </p>
                  <p className="text-sm text-white">
                    {campaignName || "Untitled"}
                  </p>
                </div>
                <div className="bg-gray-800/50 rounded-xl p-4">
                  <p className="text-xs text-gray-500 uppercase tracking-wider mb-1">
                    Type
                  </p>
                  <p className="text-sm text-white">{campaignType}</p>
                </div>
                <div className="bg-gray-800/50 rounded-xl p-4">
                  <p className="text-xs text-gray-500 uppercase tracking-wider mb-1">
                    Platforms
                  </p>
                  <p className="text-sm text-white">
                    {selectedPlatforms
                      .map((p) => p.charAt(0).toUpperCase() + p.slice(1))
                      .join(", ")}
                  </p>
                </div>
                <div className="bg-gray-800/50 rounded-xl p-4">
                  <p className="text-xs text-gray-500 uppercase tracking-wider mb-1">
                    Duration
                  </p>
                  <p className="text-sm text-white">{duration} days</p>
                </div>
                <div className="bg-gray-800/50 rounded-xl p-4">
                  <p className="text-xs text-gray-500 uppercase tracking-wider mb-1">
                    Tone
                  </p>
                  <p className="text-sm text-white">{tone}</p>
                </div>
                <div className="bg-gray-800/50 rounded-xl p-4">
                  <p className="text-xs text-gray-500 uppercase tracking-wider mb-1">
                    Hooks Selected
                  </p>
                  <p className="text-sm text-white">
                    {selectedHooks.length} hooks
                  </p>
                </div>
              </div>

              <div className="bg-gray-800/50 rounded-xl p-4">
                <p className="text-xs text-gray-500 uppercase tracking-wider mb-2">
                  Content Pillars
                </p>
                <div className="flex flex-wrap gap-2">
                  {pillars.map((p) => (
                    <span
                      key={p}
                      className="px-3 py-1 bg-purple-500/10 text-purple-400 rounded-full text-sm"
                    >
                      {p}
                    </span>
                  ))}
                </div>
              </div>

              {generating && (
                <div className="space-y-3">
                  {generationSteps.map((gs, i) => (
                    <div
                      key={gs}
                      className={`flex items-center gap-3 text-sm ${
                        i < genStep
                          ? "text-green-400"
                          : i === genStep
                          ? "text-purple-400"
                          : "text-gray-600"
                      }`}
                    >
                      {i < genStep ? (
                        <Check className="w-4 h-4" />
                      ) : i === genStep ? (
                        <Loader2 className="w-4 h-4 animate-spin" />
                      ) : (
                        <div className="w-4 h-4 rounded-full border border-gray-700" />
                      )}
                      {gs}
                    </div>
                  ))}
                  <div className="w-full bg-gray-800 rounded-full h-2 mt-2">
                    <div
                      className="bg-purple-500 h-2 rounded-full transition-all duration-500"
                      style={{
                        width: `${((genStep + 1) / generationSteps.length) * 100}%`,
                      }}
                    />
                  </div>
                </div>
              )}

              {!generating && (
                <button
                  onClick={generateCampaign}
                  className="w-full py-3 px-6 bg-purple-600 hover:bg-purple-700 text-white font-medium rounded-xl transition-colors flex items-center justify-center gap-2"
                >
                  <Sparkles className="w-4 h-4" />
                  Generate Campaign
                </button>
              )}
            </div>
          ) : (
            <div className="bg-gray-900 border border-green-500/30 rounded-2xl p-8 text-center space-y-4">
              <div className="w-16 h-16 bg-green-500/10 rounded-full flex items-center justify-center mx-auto">
                <Check className="w-8 h-8 text-green-400" />
              </div>
              <h2 className="text-2xl font-bold text-white">
                Campaign Generated!
              </h2>
              <p className="text-gray-400 max-w-md mx-auto">
                Your campaign &ldquo;{campaignName || "Untitled"}&rdquo; has
                been created with {Math.floor(duration * 0.8)} content pieces
                across {selectedPlatforms.length} platform
                {selectedPlatforms.length > 1 ? "s" : ""}.
              </p>
              <div className="flex justify-center gap-4 pt-2">
                <Link
                  href="/calendar"
                  className="flex items-center gap-2 px-5 py-2.5 bg-purple-600 hover:bg-purple-700 text-white rounded-xl text-sm font-medium transition-colors"
                >
                  <Calendar className="w-4 h-4" />
                  View Calendar
                </Link>
                <Link
                  href="/content"
                  className="flex items-center gap-2 px-5 py-2.5 bg-gray-800 hover:bg-gray-700 text-white rounded-xl text-sm font-medium transition-colors border border-gray-700"
                >
                  <Library className="w-4 h-4" />
                  Content Library
                </Link>
              </div>
            </div>
          )}
        </div>
      )}

      {/* Navigation Buttons */}
      {!(step === 3 && generated) && (
        <div className="flex justify-between mt-8">
          <button
            onClick={() => setStep(step - 1)}
            disabled={step === 0}
            className="flex items-center gap-2 px-5 py-2.5 bg-gray-800 hover:bg-gray-700 disabled:opacity-30 disabled:hover:bg-gray-800 text-white rounded-xl text-sm font-medium transition-colors border border-gray-700"
          >
            <ArrowLeft className="w-4 h-4" />
            Back
          </button>
          {step < 3 && (
            <button
              onClick={() => setStep(step + 1)}
              className="flex items-center gap-2 px-5 py-2.5 bg-purple-600 hover:bg-purple-700 text-white rounded-xl text-sm font-medium transition-colors"
            >
              Continue
              <ArrowRight className="w-4 h-4" />
            </button>
          )}
        </div>
      )}
    </div>
  );
}
