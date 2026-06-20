# JARVIS for Biology — Sprint 2 Video Script

**Target length:** ~3 minutes (~540 words spoken)
**Format:** Talking head → architecture diagram → Pareto figure → hypervolume plot → talking head
**Title (working):** Building a world model to scale our ability to explore and optimize therapeutic design space

---

## Section 1 — The setup
### We have the target. How do we design the therapeutic?
**[TALKING HEAD · 0:00–0:35]**

Two weeks ago I showed how an AI system could generate six drug targets for age-related macular degeneration in four minutes. One of them was complement gene C9 — with a hypothesis that a genetic variant was driving excess complement activation in the back of the eye, ultimately destroying photoreceptors.

This sprint I asked a different question. We have the target. How do we design the therapeutic?

---

## Section 2 — The problem
### The design space exceeds 10¹¹. Directed evolution samples it blind.
**[TALKING HEAD · 0:35–1:20]**

The therapeutic approach is gene therapy — deliver a gene directly to retinal cells that produces a protein to neutralize the complement cascade locally. Potentially a one-time treatment.

The delivery vehicle is an engineered virus called AAV2. And this is where the hard problem lives.

AAV2 has a protein shell — a capsid — of 735 amino acids. About 60 of those residues control two things simultaneously: whether the virus gets into the right cells in the eye, and whether the patient's immune system recognizes and destroys it before it ever gets there. AAV2 seroprevalence runs 50 to 70 percent in adults. You have to engineer your way out.

And the design space for that engineering problem exceeds ten-to-the-eleven possibilities. No directed evolution campaign — the standard approach in the field — will ever meaningfully cover that space. It works; we have gene therapies in the clinic. But it's blind sampling. It has no model of the space, no way to navigate tradeoffs, no principled way to enforce a safety constraint during selection.

A world model changes that.

---

## Section 3 — How it works
### Four layers, one closed loop
**[CUT TO: World model architecture diagram · 1:20–1:55]**

This is the architecture I built.

Layer one is generative — it produces about 80 capsid variants by making substitutions at the variable regions that drive tropism and antibody recognition, and by inserting peptide sequences at the site used by AAV.7m8, the best intravitreal capsid in the literature.

Layer two is feature engineering — ESM3, EvolutionaryScale's protein language model, embeds each variant sequence into a 2560-dimensional vector that captures something about its structure and likely function.

Layer three is the world model itself — a multi-output Gaussian process that learns, from every experiment, how capsid sequence relates to three outputs: transduction efficiency, antibody escape, and inflammatory response. Critically, it gives uncertainty estimates — not just predictions, but confidence.

And layer four is the closed loop — an RL policy that uses those predictions and uncertainties to pick the next candidates more intelligently than random. After each cycle, the observations feed back in, the GP refits, and the policy updates. That feedback arrow is the whole point.

---

## Section 4 — The result
### Optimizing two competing objectives under a safety constraint
**[CUT TO: Pareto frontier figure · 1:55–2:25]**

This is what it produced.

X-axis is RPE transduction efficiency — how well the capsid gets into the target cells. Y-axis is neutralizing antibody escape — what fraction of patient antibodies this capsid evades. The red shaded zone is the inflammation constraint — candidates that land there don't get selected.

The red diamond is AAV2 wildtype — the baseline. Low transduction, low escape. The green diamond is AAV.7m8 — the state of the art from directed evolution. The blue dashed line is the Pareto frontier the system found — the non-dominated candidates where you can't improve one objective without sacrificing the other.

The RL policy finds that frontier faster than random selection. By cycle four or five the gap is visible, and it keeps widening.

---

## Section 5 — Why it works
### RL policy beats random
**[CUT TO: RL-vs-random hypervolume plot · 2:25–2:40]**

And here's why the closed loop matters. This is Pareto hypervolume over cycles — a single number for how much of the optimal tradeoff space the system has covered. The orange line is the RL policy. The gray is random selection, same budget. Each cycle, observations feed back, the model refits, the policy updates — and the gap opens early and keeps widening.

---

## Section 6 — The close
### Same architecture. One leg further.
**[TALKING HEAD · 2:40–3:00]**

The wet-lab simulator is the only synthetic component in this pipeline — it's calibrated to published intravitreal AAV literature, but it's not a wet lab. Everything else is real. And when real wet-lab results are available, you swap them in. The policy interface doesn't change.

Last sprint, AI explored hypothesis space. This sprint, AI explored therapeutic design space. Same architecture. One leg further.

Full write-up, figures, and code are in the repo — link in the comments.

---

## Production notes

- **Total spoken length:** ~540 words → lands at 2:50–3:00 at a comfortable pace. The deck has a dedicated hypervolume slide, so the script now includes a short "why it works" beat for it.
- **Section ↔ slide map:** the headings above match your deck slide-for-slide — setup (slide 2), problem (slide 3), how it works (slide 4), result (slide 5), why it works (slide 6), close (slide 7).
- **Recording order:** Record all talking-head sections first against a plain background, then record the screen cuts separately while narrating live over each figure. Sync in any basic editor (CapCut, DaVinci, iMovie).
- **Visual assets needed:**
  - World model architecture diagram (`outputs/world_model.png`) → section 3
  - Pareto frontier figure (`outputs/pareto_frontier.png`) → section 4
  - RL-vs-random hypervolume plot (`outputs/rl_vs_random.png`) → section 5
- **On-screen captions:** Add burned-in captions — most LinkedIn video is watched muted. Key terms worth captioning: *C9*, *AAV2 capsid*, *10¹¹ design space*, *world model*, *Pareto frontier*, *hypervolume*.
- **Diagram timing:** The architecture diagram is dense. Consider a subtle highlight/zoom on each layer as you narrate it (purple → green → blue → orange), so the viewer's attention tracks your words.
- **Two vs. three framing:** The result headline now reads "two competing objectives under a safety constraint" to match the figure subtitle. If you ad-lib, keep inflammation as a *constraint*, not a third objective — it's a hard filter, not a Pareto axis.
