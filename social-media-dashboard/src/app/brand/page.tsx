const brandData = {
  name: "BrewCraft",
  tagline: "Artisan Coffee, Crafted with Purpose",
  colors: [
    { name: "Espresso", hex: "#3C1518" },
    { name: "Caramel", hex: "#C97D38" },
    { name: "Cream", hex: "#F5E6CC" },
    { name: "Forest", hex: "#2D4A3E" },
    { name: "Charcoal", hex: "#2B2D2F" },
  ],
  typography: {
    heading: "Playfair Display",
    body: "Inter",
  },
  voiceTraits: [
    "Warm",
    "Authentic",
    "Knowledgeable",
    "Approachable",
    "Sustainable-minded",
    "Community-driven",
    "Artisanal",
    "Inviting",
  ],
  targetAudience: {
    demographics: [
      "Age: 25-45",
      "Urban professionals",
      "Income: $50K-$120K",
      "College educated",
    ],
    psychographics: [
      "Values quality over convenience",
      "Environmentally conscious",
      "Seeks authentic experiences",
      "Active on Instagram & TikTok",
      "Willing to pay premium for craft products",
    ],
  },
  contentPillars: [
    {
      title: "Behind the Beans",
      description:
        "Origin stories, sourcing journeys, and farmer partnerships",
      percentage: 30,
    },
    {
      title: "Brew Mastery",
      description:
        "Brewing tutorials, recipes, tips, and equipment reviews",
      percentage: 25,
    },
    {
      title: "Community & Culture",
      description:
        "Customer stories, events, collaborations, and local impact",
      percentage: 25,
    },
    {
      title: "Sustainability",
      description:
        "Eco practices, packaging innovations, and impact reports",
      percentage: 20,
    },
  ],
  competitors: [
    {
      name: "Blue Bottle",
      followers: "485K",
      engagement: "3.2%",
      strength: "Premium aesthetic",
    },
    {
      name: "Stumptown",
      followers: "320K",
      engagement: "2.8%",
      strength: "Origin storytelling",
    },
    {
      name: "Counter Culture",
      followers: "198K",
      engagement: "4.1%",
      strength: "Education content",
    },
    {
      name: "Intelligentsia",
      followers: "275K",
      engagement: "3.5%",
      strength: "Community engagement",
    },
  ],
};

export default function BrandPage() {
  return (
    <div className="p-8 max-w-6xl">
      {/* Header */}
      <div className="mb-8">
        <h1 className="text-3xl font-bold text-white">{brandData.name}</h1>
        <p className="text-lg text-gray-400 mt-1">{brandData.tagline}</p>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Color Palette */}
        <div className="bg-gray-900 border border-gray-800 rounded-2xl p-6">
          <h2 className="text-lg font-semibold text-white mb-4">
            Color Palette
          </h2>
          <div className="flex gap-4">
            {brandData.colors.map((color) => (
              <div key={color.hex} className="flex flex-col items-center gap-2">
                <div
                  className="w-14 h-14 rounded-full border-2 border-gray-700"
                  style={{ backgroundColor: color.hex }}
                />
                <p className="text-xs text-white font-medium">{color.name}</p>
                <p className="text-xs text-gray-500">{color.hex}</p>
              </div>
            ))}
          </div>
        </div>

        {/* Typography */}
        <div className="bg-gray-900 border border-gray-800 rounded-2xl p-6">
          <h2 className="text-lg font-semibold text-white mb-4">Typography</h2>
          <div className="space-y-4">
            <div>
              <p className="text-xs text-gray-500 uppercase tracking-wider mb-1">
                Heading Font
              </p>
              <p className="text-2xl text-white">
                {brandData.typography.heading}
              </p>
            </div>
            <div>
              <p className="text-xs text-gray-500 uppercase tracking-wider mb-1">
                Body Font
              </p>
              <p className="text-lg text-gray-300">
                {brandData.typography.body}
              </p>
            </div>
          </div>
        </div>

        {/* Voice & Tone */}
        <div className="bg-gray-900 border border-gray-800 rounded-2xl p-6">
          <h2 className="text-lg font-semibold text-white mb-4">
            Voice &amp; Tone
          </h2>
          <div className="flex flex-wrap gap-2">
            {brandData.voiceTraits.map((trait) => (
              <span
                key={trait}
                className="px-3 py-1.5 bg-purple-500/10 text-purple-400 rounded-full text-sm font-medium"
              >
                {trait}
              </span>
            ))}
          </div>
        </div>

        {/* Target Audience */}
        <div className="bg-gray-900 border border-gray-800 rounded-2xl p-6">
          <h2 className="text-lg font-semibold text-white mb-4">
            Target Audience
          </h2>
          <div className="space-y-4">
            <div>
              <p className="text-xs text-gray-500 uppercase tracking-wider mb-2">
                Demographics
              </p>
              <div className="space-y-1">
                {brandData.targetAudience.demographics.map((item) => (
                  <p key={item} className="text-sm text-gray-300">
                    {item}
                  </p>
                ))}
              </div>
            </div>
            <div>
              <p className="text-xs text-gray-500 uppercase tracking-wider mb-2">
                Psychographics
              </p>
              <div className="space-y-1">
                {brandData.targetAudience.psychographics.map((item) => (
                  <p key={item} className="text-sm text-gray-300">
                    {item}
                  </p>
                ))}
              </div>
            </div>
          </div>
        </div>

        {/* Content Pillars */}
        <div className="lg:col-span-2 bg-gray-900 border border-gray-800 rounded-2xl p-6">
          <h2 className="text-lg font-semibold text-white mb-4">
            Content Pillars
          </h2>
          <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-4">
            {brandData.contentPillars.map((pillar) => (
              <div
                key={pillar.title}
                className="bg-gray-800/50 rounded-xl p-4"
              >
                <div className="flex items-center justify-between mb-2">
                  <p className="text-sm font-semibold text-white">
                    {pillar.title}
                  </p>
                  <span className="text-xs text-purple-400 font-bold">
                    {pillar.percentage}%
                  </span>
                </div>
                <p className="text-xs text-gray-400">{pillar.description}</p>
                <div className="mt-3 w-full bg-gray-700 rounded-full h-1.5">
                  <div
                    className="bg-purple-500 h-1.5 rounded-full"
                    style={{ width: `${pillar.percentage}%` }}
                  />
                </div>
              </div>
            ))}
          </div>
        </div>

        {/* Competitor Analysis */}
        <div className="lg:col-span-2 bg-gray-900 border border-gray-800 rounded-2xl p-6">
          <h2 className="text-lg font-semibold text-white mb-4">
            Competitor Analysis
          </h2>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="text-gray-500 text-xs uppercase tracking-wider">
                  <th className="text-left py-3 px-4">Brand</th>
                  <th className="text-left py-3 px-4">Followers</th>
                  <th className="text-left py-3 px-4">Engagement</th>
                  <th className="text-left py-3 px-4">Key Strength</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-gray-800">
                {brandData.competitors.map((comp) => (
                  <tr key={comp.name} className="hover:bg-gray-800/50">
                    <td className="py-3 px-4 text-white font-medium">
                      {comp.name}
                    </td>
                    <td className="py-3 px-4 text-gray-300">
                      {comp.followers}
                    </td>
                    <td className="py-3 px-4 text-gray-300">
                      {comp.engagement}
                    </td>
                    <td className="py-3 px-4 text-gray-400">
                      {comp.strength}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      </div>
    </div>
  );
}
