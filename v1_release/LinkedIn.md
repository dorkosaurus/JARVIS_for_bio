# Building a world model to scale our ability to explore and optimize therapeutic design space through closed-loop experimentation

Two weeks ago I demonstrated how an inference-oriented architecture could leverage agents to scale the generation and exploration of therapeutic hypothesis space.  I was able to generate six valid AMD mechanistic hypotheses in four minutes, tasks which could take a human team tens of hours.

This gave us the targets. Now we need the therapeutics.

Just like target discovery, therapeutic design space is a large combinatorial problem. And we often need to optimize against multiple objectives simultaneously with a limited number of knobs to turn.

Consider gene therapy for AMD.

AAV2 — a viral vector used in approved gene therapies — has a capsid of 735 amino acids. Roughly 60 of those residues drive two things simultaneously: site specificity (getting into the eye) and antibody recognition (neutralizing antibodies prevent the virus from ever reaching the eye).

How can we mutate those 60 residues one at a time or in combination to:

* Optimize infiltration into the retinal pigment epithelium?
* Evade antibody detection?
* Avoid triggering a systemic inflammatory response?

To do this rigorously, we'd need to explore a space of 10¹¹ possibilities. 

No directed evolution campaign will ever meaningfully cover that space. Standard campaigns generate a library, expose it to the assay, sequence what survives, repeat.  This is blind sampling- and it works. We do have gene therapies for AMD in the clinic, some of them efficacious. 

But blind sampling has no model of the space. It doesn't learn the shape of the tradeoffs, can't locate the Pareto frontier, and has no principled way to enforce a safety constraint during selection.

A world model based architecture changes the problem:

* Instead of sampling blindly, you build a probabilistic model that learns — from every experiment — how capsid sequence relates to each objective. 
* An RL policy uses that model to navigate the space rather than sample it.

That's what I built over the last two weeks.

ESM3 embeddings of real AAV2 VP1 sequences feed a multi-output Gaussian process world model. 

An RL policy — pretrained on a biologically calibrated wet-lab simulator — explores the Pareto frontier of all three objectives simultaneously. 10 cycles, 5 candidates per cycle.  The policy outperforms random selection by cycle 4–5 without a single inflammatory constraint violation on selected candidates.

The wet-lab simulator is the only synthetic component in this architecture and we can swap in real wet lab results when we have them.  The policy interface won't have to change. 

Full write-up, figures, and reproducible code in the comments.

What's the hardest design tradeoff you've navigated in your area of biology?

#AIxScience #drugdiscovery #generativebiology #BuildingJARVIS
