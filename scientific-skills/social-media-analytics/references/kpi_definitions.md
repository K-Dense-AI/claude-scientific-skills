# Social Media KPI Definitions

This reference provides comprehensive definitions, formulas, benchmarks, and improvement strategies for every KPI tracked by the social-media-analytics skill.

---

## Awareness Metrics

### Reach

- **Definition**: The number of unique users who saw your content at least once during the reporting period.
- **Formula**: `Reach = count of unique users served the content`
- **What It Measures**: The breadth of your content distribution. High reach means your content is surfacing to a wide audience.
- **Ranges**:
  - Good: >15% of followers per post (organic)
  - Average: 5-15% of followers per post
  - Poor: <5% of followers per post
- **How to Improve**: Post at peak activity hours, use relevant hashtags, encourage shares, experiment with Reels/Shorts which receive algorithmic reach boosts.
- **Platform Nuances**:
  - **Instagram**: Reels reach is typically 2-3x higher than static posts. Explore page placement dramatically increases reach.
  - **Facebook**: Page reach has declined to 5-7% organically. Groups and Reels extend reach.
  - **LinkedIn**: Organic reach is strong (8-12% of connections) especially for text-only and document posts.
  - **TikTok**: Reach is decoupled from follower count; the For You Page can push reach to millions for any account.

### Impressions

- **Definition**: The total number of times your content was displayed, including multiple views by the same user.
- **Formula**: `Impressions = total display count (includes repeat views)`
- **What It Measures**: Overall content visibility. Impressions are always >= reach. A high impression-to-reach ratio means people are seeing your content more than once.
- **Ranges**:
  - Good: Impressions/Reach ratio of 1.5-3.0 (indicates repeat exposure)
  - Average: Ratio of 1.1-1.5
  - Poor: Ratio near 1.0 (each person sees it only once, weak retention)
- **How to Improve**: Create content worth revisiting (infographics, reference material). Use carousel formats that users swipe through multiple times. Pin top posts.
- **Platform Nuances**:
  - **Instagram**: Stories generate high impressions from repeat viewers. Carousels inflate impressions as each slide counts.
  - **Facebook**: Shared posts in groups generate secondary impression spikes.
  - **LinkedIn**: Re-engagement from comments on older posts can sustain impressions over several days.
  - **TikTok**: Loop-friendly short videos can generate impression-to-reach ratios above 3.0.

### Share of Voice (SOV)

- **Definition**: Your brand's share of total mentions or conversations within your industry or competitive set.
- **Formula**: `SOV = (brand mentions / total industry mentions) * 100`
- **What It Measures**: Market visibility relative to competitors. Indicates how dominant your brand is in social conversations.
- **Ranges**:
  - Good: >25% in your category
  - Average: 10-25%
  - Poor: <10%
- **How to Improve**: Increase posting volume, run hashtag campaigns, partner with influencers, engage in trending conversations relevant to your industry.
- **Platform Nuances**:
  - **Twitter/X**: Most relevant platform for SOV due to public conversation volume.
  - **LinkedIn**: SOV matters for B2B thought leadership.
  - **Instagram/TikTok**: Harder to measure; use branded hashtag volume as a proxy.

### Brand Mention Rate

- **Definition**: The frequency at which your brand is mentioned by others (not your own posts) across social platforms.
- **Formula**: `Brand Mention Rate = brand mentions / total posts in category`
- **What It Measures**: Organic buzz and word-of-mouth traction.
- **Ranges**:
  - Good: Consistent week-over-week growth in mentions
  - Average: Stable mention volume
  - Poor: Declining mentions or mentions only from owned accounts
- **How to Improve**: Create shareable content, run UGC campaigns, engage in community conversations, launch referral incentives.
- **Platform Nuances**:
  - Mentions are easier to track on Twitter/X and LinkedIn (public APIs). Instagram mentions require monitoring tags and stories. TikTok mentions are primarily via duets, stitches, and comments.

---

## Engagement Metrics

### Engagement Rate

- **Definition**: The percentage of people who interacted with your content (likes, comments, shares, saves, clicks) relative to the audience size.
- **Formula (by reach)**: `ER = (total engagements / reach) * 100`
- **Formula (by followers)**: `ER = (total engagements / followers) * 100`
- **What It Measures**: How compelling your content is to the audience that sees it. The single most important metric for content quality.
- **Ranges** (by followers):
  - Good: >3.5%
  - Average: 1.0-3.5%
  - Poor: <1.0%
- **How to Improve**: Ask questions in captions, use strong CTAs, post content that triggers emotional responses (humor, inspiration, controversy), respond to every comment within the first hour.
- **Platform Nuances**:
  - **Instagram**: Average ER is 1.5-3.0% (higher for <10K follower accounts). Reels typically achieve 1.5-2x the ER of static posts.
  - **Facebook**: Average ER is 0.5-1.0%. Video and live content outperform links and text.
  - **LinkedIn**: Average ER is 2.0-4.0%. Document/carousel posts achieve the highest ER (3-5%).
  - **TikTok**: Average ER is 4.0-8.0%. The highest of all major platforms due to algorithmic distribution.

### Comment Rate

- **Definition**: The number of comments per post relative to followers or reach.
- **Formula**: `Comment Rate = (comments / followers) * 100`
- **What It Measures**: Depth of engagement. Comments require more effort than likes, indicating stronger audience connection.
- **Ranges**:
  - Good: >0.5% of followers
  - Average: 0.1-0.5%
  - Poor: <0.1%
- **How to Improve**: End posts with open-ended questions, create debate-worthy content, reply to comments to spark threads, use "caption this" or "fill in the blank" formats.
- **Platform Nuances**:
  - **Instagram**: Comment pods can artificially inflate; watch for genuine vs. manufactured comments.
  - **LinkedIn**: Long-form posts and personal stories drive the highest comment rates.
  - **TikTok**: "Reply to comment" videos create engagement loops that boost overall comment rates.

### Comment Sentiment

- **Definition**: The distribution of positive, neutral, and negative sentiment across comments on your content.
- **Formula**: `Positive Ratio = positive comments / total comments` (same for neutral and negative)
- **What It Measures**: Audience perception and brand health. A high negative ratio can indicate a PR issue even if total engagement is strong.
- **Ranges**:
  - Good: >70% positive, <10% negative
  - Average: 50-70% positive, 10-20% negative
  - Poor: <50% positive or >20% negative
- **How to Improve**: Address negative feedback promptly, moderate toxic comments, create content that aligns with audience values, avoid controversial topics that don't align with brand positioning.
- **Platform Nuances**:
  - **TikTok**: Comment sections tend to be more casual and humorous; "negative" comments may be ironic.
  - **LinkedIn**: Comments are generally more professional; negative sentiment often signals genuine disagreement.
  - **Instagram**: Emoji-heavy responses require sentiment models that handle non-text input.

### Share/Repost Rate

- **Definition**: The number of shares, reposts, or retweets per post relative to followers.
- **Formula**: `Share Rate = (shares / followers) * 100`
- **What It Measures**: Content virality potential. Shares extend reach beyond your existing audience.
- **Ranges**:
  - Good: >1.0%
  - Average: 0.3-1.0%
  - Poor: <0.3%
- **How to Improve**: Create highly informative or emotionally resonant content. Infographics, data visualizations, and "tag someone who..." posts drive shares. Make content easy to share (clean visuals, clear value proposition).
- **Platform Nuances**:
  - **Facebook**: Shares are the primary distribution mechanism for organic reach.
  - **LinkedIn**: Reposts plus "shared with thoughts" posts both count. Industry data and contrarian takes drive shares.
  - **TikTok**: Shares include DM sends (private) and reposts (public). Humor and "relatable" content dominates shares.
  - **Instagram**: Shares to Stories and DMs are the main vectors. Instagram does not show share counts publicly.

### Save/Bookmark Rate

- **Definition**: The number of saves or bookmarks per post relative to reach.
- **Formula**: `Save Rate = (saves / reach) * 100`
- **What It Measures**: Content utility. Users save content they want to return to, signaling high-value material. Saves are a strong positive signal to recommendation algorithms.
- **Ranges**:
  - Good: >3.0% of reach
  - Average: 1.0-3.0%
  - Poor: <1.0%
- **How to Improve**: Create reference-worthy content (checklists, templates, how-to guides, data tables). Use carousels with step-by-step instructions. Include "save this for later" CTAs.
- **Platform Nuances**:
  - **Instagram**: Saves are weighted heavily in the algorithm. Carousels and infographics are the most saved formats.
  - **TikTok**: Bookmarks signal educational value; tutorials and recipes are highly saved.
  - **LinkedIn**: Saves are less visible but exist; long-form how-to posts are saved most.

### Amplification Rate

- **Definition**: The ratio of shares to followers, measuring how much your audience amplifies your message to their own networks.
- **Formula**: `Amplification Rate = (shares / followers) * 100`
- **What It Measures**: How effectively your content turns followers into distributors.
- **Ranges**:
  - Good: >2.0%
  - Average: 0.5-2.0%
  - Poor: <0.5%
- **How to Improve**: Create content that makes sharers look good to their own audience (insights, humor, breaking news). Align with trending topics.
- **Platform Nuances**: Similar to Share Rate. Most relevant on Twitter/X and LinkedIn where sharing is frictionless.

### Virality Rate

- **Definition**: The ratio of shares to impressions, measuring how likely someone who sees your content is to share it.
- **Formula**: `Virality Rate = (shares / impressions) * 100`
- **What It Measures**: Per-impression share probability. A purer measure of content shareability than amplification rate.
- **Ranges**:
  - Good: >1.0%
  - Average: 0.3-1.0%
  - Poor: <0.3%
- **How to Improve**: Optimize the first 3 seconds of video for immediate hook. Use surprising statistics or counterintuitive claims. Create content with a strong "I need to show this to someone" impulse.
- **Platform Nuances**:
  - **TikTok**: Virality rate benchmarks are higher due to the platform's share-centric culture.
  - **LinkedIn**: Virality is lower in absolute terms but high-quality B2B content can achieve spikes.

---

## Content Metrics

### Click-Through Rate (CTR)

- **Definition**: The percentage of impressions that resulted in a link click.
- **Formula**: `CTR = (link clicks / impressions) * 100`
- **What It Measures**: How effectively your content drives traffic to an external destination (website, landing page, product page).
- **Ranges**:
  - Good: >2.0%
  - Average: 0.8-2.0%
  - Poor: <0.8%
- **How to Improve**: Use compelling CTAs, place links in the first comment (LinkedIn), use link stickers in Stories, write benefit-driven copy that creates curiosity gaps.
- **Platform Nuances**:
  - **Facebook**: Link posts have lower organic reach than native content. Use video or image posts with links in comments.
  - **Instagram**: Links are only clickable in bio, Stories (sticker), and DMs. Linktree usage is common.
  - **LinkedIn**: Links in post body reduce reach by ~40%. Place links in the first comment instead.
  - **TikTok**: Links are only available via bio or TikTok Shop. CTR is measured differently.

### Video View Rate

- **Definition**: The percentage of impressions that resulted in a video view (platform-defined threshold).
- **Formula**: `Video View Rate = (video views / impressions) * 100`
- **What It Measures**: How effectively your video thumbnail and first frame capture attention.
- **Ranges**:
  - Good: >30%
  - Average: 15-30%
  - Poor: <15%
- **How to Improve**: Use attention-grabbing thumbnails, start with a hook in the first 1-2 seconds, add captions for sound-off viewing.
- **Platform Nuances**:
  - **Facebook**: A "view" counts at 3 seconds. Autoplay inflates view counts.
  - **Instagram**: Reels count a view at 3 seconds (or a loop for short videos).
  - **LinkedIn**: A view counts at 2 seconds with at least 50% of the video on screen.
  - **TikTok**: A view is counted as soon as the video starts playing.

### Average Watch Duration

- **Definition**: The mean number of seconds viewers watch your video content.
- **Formula**: `Avg Watch Duration = total watch time / total video views`
- **What It Measures**: Content retention and quality. Algorithms heavily weight watch duration for distribution decisions.
- **Ranges**:
  - Good: >50% of video length
  - Average: 25-50%
  - Poor: <25%
- **How to Improve**: Front-load value, use pattern interrupts every 3-5 seconds, create narrative tension, keep videos as short as the content allows.
- **Platform Nuances**:
  - **TikTok**: Watch duration is the single most important algorithm signal. Aim for >80% on short-form (<30s) videos.
  - **Instagram Reels**: Similar to TikTok; completion rate and re-watches are key.
  - **LinkedIn**: Longer watch times are acceptable (60-90s is the sweet spot) for educational content.

### Video Completion Rate

- **Definition**: The percentage of viewers who watched the video to the end.
- **Formula**: `Completion Rate = (viewers who watched to end / total viewers) * 100`
- **What It Measures**: End-to-end engagement. High completion rates signal compelling content and trigger algorithmic amplification.
- **Ranges**:
  - Good: >40%
  - Average: 20-40%
  - Poor: <20%
- **How to Improve**: Keep videos concise, promise a payoff at the end ("wait for it"), use countdown or list formats, tease the conclusion in the opening.
- **Platform Nuances**:
  - **TikTok**: Short videos (<15s) routinely achieve 60-80% completion. Longer videos need stronger hooks.
  - **Instagram Reels**: Completion rates above 50% significantly boost Explore page distribution.
  - **Facebook**: Completion rates are lower overall due to passive scrolling behavior.

### Story Completion Rate

- **Definition**: The percentage of viewers who watched all frames of your Story (from first to last).
- **Formula**: `Story Completion Rate = (viewers of last frame / viewers of first frame) * 100`
- **What It Measures**: How engaging your multi-frame Story content is. Drop-off between frames indicates where interest wanes.
- **Ranges**:
  - Good: >70%
  - Average: 50-70%
  - Poor: <50%
- **How to Improve**: Keep Stories to 3-7 frames, use interactive stickers (polls, questions), maintain visual consistency, place the most engaging content in the first 2 frames.
- **Platform Nuances**:
  - **Instagram**: Most Stories see a 15-25% drop-off after the first frame. Interactive stickers increase completion by 10-15%.
  - **Facebook**: Story viewership is lower but completion rates are comparable to Instagram.

### Story Reply Rate

- **Definition**: The number of direct replies to your Story relative to Story views.
- **Formula**: `Story Reply Rate = (replies / story views) * 100`
- **What It Measures**: Conversational engagement. Story replies are private and indicate strong personal connection with the content.
- **Ranges**:
  - Good: >2.0%
  - Average: 0.5-2.0%
  - Poor: <0.5%
- **How to Improve**: Ask direct questions, use the question sticker, create "fill in the blank" prompts, share behind-the-scenes content that invites personal responses.
- **Platform Nuances**:
  - **Instagram**: Replies are DMs; high reply rates build 1:1 relationships. Algorithms may boost accounts with high reply rates.

---

## Growth Metrics

### Follower Growth Rate

- **Definition**: The net change in followers over a period, expressed as a percentage of starting followers.
- **Formula**: `FGR = ((new followers - unfollows) / starting followers) * 100`
- **What It Measures**: How quickly your audience is growing (or shrinking). Combines acquisition and retention.
- **Ranges** (monthly):
  - Good: >3.0%
  - Average: 1.0-3.0%
  - Poor: <1.0% or negative
- **How to Improve**: Maintain consistent posting schedule, cross-promote across platforms, collaborate with complementary accounts, run giveaways, optimize bio and profile for conversions.
- **Platform Nuances**:
  - **TikTok**: Growth can be explosive (10-50% monthly) due to viral potential but also volatile.
  - **LinkedIn**: Slower growth (1-3%) but followers are typically higher-value for B2B.
  - **Instagram**: Growth has slowed for most accounts; Reels are the primary growth driver.

### Follower-to-Engagement Ratio

- **Definition**: Engagement rate contextualized by follower count tier.
- **Formula**: `F2E Ratio = engagement rate / expected engagement rate for follower tier`
- **What It Measures**: Whether engagement is proportional to audience size. Smaller accounts naturally have higher ER; this metric normalizes for that.
- **Ranges**:
  - Good: >1.2 (outperforming your tier)
  - Average: 0.8-1.2 (on par)
  - Poor: <0.8 (underperforming your tier)
- **How to Improve**: Focus on audience quality over quantity, prune inactive followers if possible, create niche content that resonates deeply with your specific audience.
- **Platform Nuances**: Apply tier-specific baselines: <10K followers expect 3-5% ER; 10K-100K expect 1.5-3%; 100K-1M expect 0.8-1.5%; >1M expect 0.5-1.0%.

### Profile Visit Rate

- **Definition**: The number of profile visits per impression.
- **Formula**: `Profile Visit Rate = (profile visits / impressions) * 100`
- **What It Measures**: How effectively your content piques curiosity about your brand or account.
- **Ranges**:
  - Good: >5.0%
  - Average: 2.0-5.0%
  - Poor: <2.0%
- **How to Improve**: Use branded content that sparks curiosity, include "link in bio" CTAs, ensure your username is visible in shared content.
- **Platform Nuances**:
  - **Instagram**: Profile visits are a key top-of-funnel metric. Reels from the Explore page drive the most profile visits.
  - **TikTok**: Profile visits spike after viral videos; ensure your bio and pinned content convert visitors.

### Website Click Rate

- **Definition**: The number of website clicks from your profile relative to profile visits.
- **Formula**: `Website Click Rate = (website clicks / profile visits) * 100`
- **What It Measures**: How effectively your profile converts curious visitors into website traffic.
- **Ranges**:
  - Good: >5.0%
  - Average: 2.0-5.0%
  - Poor: <2.0%
- **How to Improve**: Keep your bio link updated and relevant, use Linktree or similar for multiple destinations, add a clear value proposition near the link, use "link in bio" CTAs in posts.
- **Platform Nuances**:
  - **Instagram**: Single link in bio (or Linktree). Story link stickers bypass the profile visit step.
  - **LinkedIn**: Website clicks from profile are tracked; also consider link clicks from posts.
  - **TikTok**: TikTok Shop links compete with bio links for click attention.

---

## Conversion Metrics

### Conversion Rate

- **Definition**: The percentage of link clicks that result in a desired action (purchase, signup, download).
- **Formula**: `Conversion Rate = (conversions / link clicks) * 100`
- **What It Measures**: The quality of traffic from social media and the effectiveness of landing pages.
- **Ranges**:
  - Good: >3.0%
  - Average: 1.0-3.0%
  - Poor: <1.0%
- **How to Improve**: Ensure message match between social content and landing page, reduce friction in the conversion flow, use social proof, retarget engaged users.
- **Platform Nuances**:
  - **Facebook/Instagram**: Conversion tracking via Meta Pixel. In-app checkout increases conversion rates.
  - **LinkedIn**: B2B conversion rates are lower (0.5-2%) but deal values are higher.
  - **TikTok**: TikTok Shop integration dramatically improves conversion rates for e-commerce.

### Cost Per Click (CPC)

- **Definition**: The average cost of each click on a paid social ad.
- **Formula**: `CPC = total ad spend / total clicks`
- **What It Measures**: Advertising efficiency in driving traffic.
- **Ranges**:
  - Good: <$0.50 (varies significantly by industry and platform)
  - Average: $0.50-$2.00
  - Poor: >$2.00
- **How to Improve**: Improve ad relevance scores, refine audience targeting, test multiple creatives, use lookalike audiences, bid strategically.
- **Platform Nuances**:
  - **Facebook**: Average CPC is $0.50-$1.50. Varies by industry (finance ~$3.50, retail ~$0.45).
  - **Instagram**: Average CPC is $0.70-$2.00. Story ads tend to have lower CPC than feed ads.
  - **LinkedIn**: Average CPC is $2.00-$5.00. Higher due to B2B targeting precision.
  - **TikTok**: Average CPC is $0.30-$1.00. Generally the lowest CPC of major platforms.

### Cost Per Engagement (CPE)

- **Definition**: The average cost of each engagement (like, comment, share, save) on paid content.
- **Formula**: `CPE = total ad spend / total engagements`
- **What It Measures**: How cost-effectively paid content generates interactions.
- **Ranges**:
  - Good: <$0.10
  - Average: $0.10-$0.50
  - Poor: >$0.50
- **How to Improve**: Use engagement-optimized ad objectives, create content that naturally encourages interaction, target warm audiences who already know your brand.
- **Platform Nuances**:
  - **Facebook**: CPE averages $0.10-$0.30.
  - **Instagram**: CPE averages $0.15-$0.40.
  - **LinkedIn**: CPE averages $0.50-$1.50 (higher due to professional audience).
  - **TikTok**: CPE averages $0.05-$0.20 (lowest due to high engagement norms).

### Return on Ad Spend (ROAS)

- **Definition**: Revenue generated for every dollar spent on social media advertising.
- **Formula**: `ROAS = revenue from ads / ad spend`
- **What It Measures**: Direct financial return of paid social campaigns.
- **Ranges**:
  - Good: >4.0x (four dollars revenue per dollar spent)
  - Average: 2.0-4.0x
  - Poor: <2.0x
- **How to Improve**: Optimize conversion funnels, use retargeting, test creative variations, focus spend on highest-performing audiences and placements, improve landing page conversion rates.
- **Platform Nuances**:
  - **Facebook/Instagram**: Mature ad platforms with strong ROAS for e-commerce (3-5x average).
  - **TikTok**: ROAS is improving rapidly; TikTok Shop integration drives 4-6x for products with viral potential.
  - **LinkedIn**: ROAS is harder to measure due to longer B2B sales cycles; use lead value estimates.

### Content ROI

- **Definition**: The total return on all social media investment (not just ads), including time, tools, and content production costs.
- **Formula**: `Content ROI = ((revenue attributed to social - total social costs) / total social costs) * 100`
- **What It Measures**: The overall business value of your social media program.
- **Ranges**:
  - Good: >200% (3x return)
  - Average: 50-200%
  - Poor: <50% or negative
- **How to Improve**: Reduce production costs by repurposing content, focus on high-ROI platforms, cut low-performing content types, invest in organic reach to reduce paid dependency.
- **Platform Nuances**: ROI varies dramatically by business model. E-commerce can track direct sales; B2B often relies on lead attribution models; brand awareness campaigns use proxy metrics (reach, SOV).
