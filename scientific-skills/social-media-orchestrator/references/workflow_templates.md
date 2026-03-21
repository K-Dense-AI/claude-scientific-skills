# Workflow Templates

End-to-end workflow templates for common social media campaign scenarios. Each template provides numbered steps with skill invocations, expected outputs, and decision points.

---

## 1. 30-Day Brand Launch Campaign

**Use when**: Launching a new brand, product, or major rebrand on social media for the first time.

**Duration**: 1 week preparation + 4 weeks active campaign

### Preparation Week (Week 0)

1. **Collect brand inputs**
   - Gather: website URL, any existing social handles, brand guidelines (if available), target audience description, campaign budget, launch date
   - Decision point: Is there an existing brand identity to work from, or starting from scratch?

2. **Run brand discovery**
   - Invoke: `brand-identity-miner`
   - Input: Website URL + social handles
   - Output: `brand_dna.json`
   - Review: Validate brand DNA with stakeholder before proceeding

3. **Generate campaign plan**
   - Run: `python scripts/orchestrate_campaign.py --brand-url <URL> --campaign-type product_launch --platforms instagram,linkedin --duration-days 30 --output-dir ./launch_campaign`
   - Output: `campaign_plan.json` with full phase breakdown

4. **Generate campaign brief**
   - Run: `python scripts/brand_brief_generator.py --brand-dna ./launch_campaign/brand_dna.json --campaign-type product_launch --output ./launch_campaign/campaign_brief.json`
   - Output: Campaign-specific brief with tone adjustments and content focus areas

5. **Create content calendar**
   - Invoke: `social-media-scheduler` for calendar generation
   - Input: campaign_plan.json + brand_dna.json
   - Output: `content_calendar.json` with 30 days of scheduled slots
   - Decision point: Approve calendar structure before content creation begins

### Week 1: Tease Phase

6. **Generate tease content briefs**
   - Invoke: `viral-content-creator`
   - Generate 5-7 tease-phase content briefs (behind-the-scenes, countdowns, sneak peeks)
   - Output: `content_briefs/tease_001.json` through `tease_007.json`

7. **Produce tease content**
   - Invoke: `social-media-copywriting` for captions
   - Invoke: `social-media-visual-generator` for visuals
   - Output: Copy + visuals for all tease posts

8. **Review and approve tease batch**
   - Present to user for approval
   - Brand consistency check against brand_dna.json
   - Decision point: Approve, revise, or reject individual posts

9. **Schedule tease content**
   - Invoke: `social-media-scheduler` publish workflow
   - Schedule posts for Week 1 dates
   - Output: publish_log entries

### Week 2: Tease Intensifies + Pre-Launch

10. **Generate pre-launch content**
    - 5-7 briefs: polls, stronger teasers, influencer seeding content
    - Content should build toward the launch moment

11. **Produce and approve pre-launch content**
    - Same flow as steps 7-8
    - Decision point: Is audience engagement sufficient? Adjust if needed.

12. **Monitor tease engagement**
    - Invoke: `social-media-analytics` for early metrics
    - Assess: Which tease posts performed best?
    - Adjust: Lean into top-performing formats for remaining content

### Week 3: Launch Phase

13. **Generate launch content**
    - Hero launch post, product demo carousel, feature highlight series
    - Launch-day Stories/Reels
    - 10-12 content pieces for heavy launch week

14. **Produce all launch assets**
    - Priority on hero visuals and demo content
    - Multiple format variations (feed, Stories, Reels)

15. **Execute launch day**
    - Publish launch content on schedule
    - Monitor engagement in real-time
    - Respond to comments and DMs within 1 hour
    - Decision point: If launch post underperforms, boost with paid promotion

### Week 4: Sustain Phase

16. **Generate sustain content**
    - Tutorials, customer testimonials, FAQ posts, UGC reshares
    - Shift from excitement to value and proof

17. **Produce and publish sustain content**
    - Same production flow
    - Incorporate any UGC that appeared during launch

18. **Run 7-day analytics check**
    - Invoke: `social-media-analytics`
    - Identify top and bottom performers
    - Adjust remaining content based on insights

19. **Final campaign report**
    - Run: `python scripts/campaign_report.py --campaign-plan ./launch_campaign/campaign_plan.json --analytics ./launch_campaign/analytics/report_30d.json --output ./launch_campaign/campaign_report.md --format md`
    - Output: Complete campaign performance report
    - Decision point: Continue with ongoing content, or plan next campaign cycle?

---

## 2. Weekly Content Production Workflow

**Use when**: Running ongoing social media content on a weekly production cycle.

**Duration**: Repeats weekly

### Monday: Planning

1. **Review previous week's analytics**
   - Invoke: `social-media-analytics`
   - Identify top 3 posts by engagement, bottom 3 posts
   - Note: What topics, formats, and times performed best?

2. **Check trending content and industry news**
   - WebSearch for trending topics in the brand's niche
   - Note any trends worth participating in this week

3. **Generate weekly content briefs**
   - Invoke: `viral-content-creator` with brand_dna.json + trend insights
   - Generate 5-7 content briefs for the week
   - Apply 60/30/10 mix: 3-4 value, 2 curated, 1 promotional
   - Output: Weekly content briefs

### Tuesday: Content Creation

4. **Produce copy for all weekly posts**
   - Invoke: `social-media-copywriting`
   - Generate captions, CTAs, hashtag sets for each post
   - Output: `captions/week_XX_*.json`

5. **Produce visuals for all weekly posts**
   - Invoke: `social-media-visual-generator`
   - Generate images, carousels, video storyboards
   - Output: `visuals/week_XX_*.png`

### Wednesday: Review

6. **Present content batch for approval**
   - Display all posts (visual + copy) for user review
   - Brand consistency check
   - Quality checklist verification

7. **Revision loop**
   - Make requested changes
   - Decision point: All posts approved? Proceed to scheduling.

### Thursday: Scheduling

8. **Schedule all approved posts**
   - Invoke: `social-media-scheduler`
   - Schedule across platforms at optimal times
   - Verify scheduling confirmation for each post

### Friday: Community Management

9. **Engage with audience**
   - Respond to comments from the week
   - Engage with community content
   - Identify UGC opportunities for next week

10. **Quick performance check**
    - Review early engagement on posts published this week
    - Note any immediate adjustments needed

---

## 3. Monthly Analytics Review Workflow

**Use when**: Conducting a comprehensive monthly performance review and optimization cycle.

**Duration**: Monthly, 1-2 hours

### Data Collection

1. **Pull 30-day analytics**
   - Invoke: `social-media-analytics` for all active platforms
   - Collect: engagement rates, reach, follower growth, top/bottom posts
   - Output: `analytics/report_30d.json`

2. **Compile content inventory**
   - List all content published in the month
   - Categorize by: platform, format, content pillar, day/time

### Analysis

3. **Identify performance patterns**
   - Top 5 posts by engagement: What do they share? (format, topic, time, visual style)
   - Bottom 5 posts: What underperformed? Why?
   - Best performing day/time combinations
   - Best performing content pillar

4. **Benchmark against goals**
   - Compare actual KPIs against campaign targets
   - Calculate: engagement rate trend, follower growth rate, reach trend
   - Decision point: Are we on track? Ahead? Behind?

5. **Competitor check**
   - Quick review of competitor activity this month
   - Note any content strategies to adapt or differentiate from

### Optimization

6. **Generate optimization recommendations**
   - Content mix adjustments (more of what works, less of what doesn't)
   - Posting time refinements based on data
   - Format emphasis changes (e.g., increase Reels if they outperform)
   - Content pillar rebalancing if needed
   - Output: `optimization_recommendations.json`

7. **Update content strategy**
   - Apply recommendations to next month's content planning
   - Update brand_dna.json if brand voice or visual style should evolve
   - Adjust posting cadence if data supports it

8. **Generate monthly report**
   - Run: `python scripts/campaign_report.py`
   - Output: Monthly performance report for stakeholders

---

## 4. Competitor Response Workflow

**Use when**: A competitor launches a major campaign, product, or initiative that requires a strategic response.

**Duration**: 1-3 days rapid response

### Day 1: Intelligence Gathering

1. **Analyze competitor campaign**
   - Invoke: `brand-identity-miner` on competitor URL/handles
   - Document: Campaign theme, messaging, platforms, content types, engagement
   - Note: What audience pain point are they addressing?

2. **Assess impact on your brand**
   - Is the competitor directly targeting your audience?
   - Are they making claims your brand should address?
   - Is there an opportunity to differentiate?
   - Decision point: React (direct response), Differentiate (indirect), or Ignore?

### Day 1-2: Strategy

3. **Define response approach**
   - **React**: Address competitor claims directly (use cautiously)
   - **Differentiate**: Highlight unique strengths without naming competitor
   - **Ride the wave**: If competitor generated industry buzz, join the conversation
   - **Counter-program**: Launch your own campaign to recapture attention

4. **Generate response content briefs**
   - Invoke: `viral-content-creator` with competitor context
   - Generate 3-5 rapid-response content pieces
   - Focus on brand strengths and unique value propositions

### Day 2-3: Execution

5. **Produce response content**
   - Invoke: `social-media-copywriting` + `social-media-visual-generator`
   - Expedited production for 3-5 posts

6. **Approve and publish**
   - Streamlined review (stakeholder fast-track approval)
   - Publish immediately or within 24 hours
   - Decision point: Boost with paid promotion for faster reach?

7. **Monitor response**
   - Track engagement on response content vs. normal posts
   - Monitor sentiment in comments
   - Adjust if response is not landing as intended

---

## 5. Crisis Communication Workflow

**Use when**: A brand faces a PR crisis, negative viral content, customer complaint escalation, or reputational threat on social media.

**Duration**: Immediate response, 1-7 days active management

### Immediate (Within 1-2 Hours)

1. **Assess the situation**
   - What happened? Scope of negative content.
   - Which platforms are affected?
   - How fast is it spreading? (volume of mentions, shares)
   - Is mainstream media involved?
   - Decision point: Is this a genuine crisis or a minor complaint?

2. **Pause scheduled content**
   - Immediately pause all scheduled posts across platforms
   - Review upcoming content for anything tone-deaf given the crisis
   - Decision point: Which scheduled posts can proceed and which must be held?

3. **Draft initial response**
   - Acknowledge the situation (do not ignore)
   - Express concern/empathy without admitting fault prematurely
   - Promise to investigate and provide updates
   - Route to appropriate stakeholders for approval
   - Tone: Empathetic, transparent, accountable

### Day 1: Active Response

4. **Publish official statement**
   - Single, consistent statement across all active platforms
   - Pin the statement to profile where possible
   - Reply to major threads/comments with link to statement

5. **Monitor conversation**
   - Track mentions, sentiment, volume
   - Identify key voices (journalists, influencers, affected customers)
   - Respond to direct questions with factual, empathetic answers
   - Do NOT engage with trolls or bad-faith actors

6. **Internal coordination**
   - Ensure customer support teams have consistent messaging
   - Prepare FAQ for common questions
   - Update stakeholders regularly

### Days 2-7: Recovery

7. **Provide updates as promised**
   - Share what actions have been taken
   - Demonstrate accountability and corrective measures
   - Decision point: Is the situation de-escalating? Adjust tone accordingly.

8. **Gradually resume regular content**
   - Start with value-first content (educational, helpful)
   - Avoid promotional content for at least 48-72 hours after crisis peaks
   - Show actions taken, not just words

9. **Post-crisis review**
   - What happened and why?
   - How effective was the response?
   - What processes need to change to prevent recurrence?
   - Update crisis communication playbook
   - Decision point: Is there lasting reputational damage that requires an ongoing repair campaign?

---

## 6. Seasonal Campaign Template

**Use when**: Planning a campaign around a holiday, season, or cultural moment.

**Duration**: 2-4 weeks of planning + 1-3 weeks active campaign

### 4 Weeks Before: Planning

1. **Select seasonal moment**
   - Identify the holiday, season, or cultural event
   - Define relevance to the brand (authentic connection required)
   - Reference: `campaign_frameworks.md` seasonal calendar

2. **Generate campaign concept**
   - Invoke: `viral-content-creator` with seasonal theme + brand_dna.json
   - Define: campaign hashtag, visual theme, key messages
   - Decision point: Is the concept authentic to the brand? Avoid forced connections.

3. **Create campaign plan**
   - Run: `python scripts/orchestrate_campaign.py --campaign-type seasonal --duration-days 14-21`
   - Define content calendar with seasonal content mix
   - Plan any promotions, offers, or limited-time content

### 2 Weeks Before: Production

4. **Generate all content briefs**
   - Invoke: `viral-content-creator`
   - 10-15 content briefs for the campaign period
   - Mix: themed value (35%), promotional (30%), emotional/storytelling (20%), interactive (15%)

5. **Produce content**
   - Invoke: `social-media-copywriting` for seasonal captions
   - Invoke: `social-media-visual-generator` for themed visuals
   - Ensure visuals adapt brand palette to seasonal theme (without losing brand identity)

6. **Review and approve full campaign**
   - Present complete content set for approval
   - Verify seasonal theme is consistent across all posts
   - Decision point: Are promotional offers finalized and approved?

### 1 Week Before: Scheduling

7. **Schedule all content**
   - Invoke: `social-media-scheduler`
   - Schedule based on seasonal posting calendar
   - Include any automated responses or Story templates

8. **Prepare real-time content templates**
   - Create templates for day-of Stories, live content
   - Prepare UGC reshare templates
   - Set up monitoring for campaign hashtag

### Campaign Active: Execution

9. **Launch campaign**
   - Publish kickoff post
   - Activate any paid promotion
   - Monitor hashtag and engagement

10. **Daily management**
    - Post scheduled content
    - Engage with audience responses
    - Share UGC and community content
    - Monitor competitor seasonal campaigns

11. **Mid-campaign check**
    - Invoke: `social-media-analytics`
    - Assess performance against goals
    - Decision point: Adjust remaining content based on performance?

### Post-Campaign: Wrap-Up

12. **Campaign wrap-up content**
    - Thank-you post, highlight reel, or recap
    - Transition back to regular content calendar

13. **Generate performance report**
    - Run: `python scripts/campaign_report.py`
    - Compare against seasonal campaign benchmarks
    - Document lessons for the next seasonal campaign

---

## Workflow Selection Guide

| Situation | Workflow | Urgency |
|-----------|----------|---------|
| New brand or product launching | 30-Day Brand Launch | Planned (4-6 weeks lead) |
| Ongoing content production | Weekly Content Production | Recurring (every week) |
| Monthly review and optimization | Monthly Analytics Review | Recurring (every month) |
| Competitor makes a big move | Competitor Response | Rapid (1-3 days) |
| PR crisis or negative virality | Crisis Communication | Immediate (hours) |
| Holiday or seasonal moment | Seasonal Campaign | Planned (4 weeks lead) |
