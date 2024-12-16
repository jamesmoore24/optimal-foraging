# Optimal Foraging in Memory Retrieval: Evaluating Random Walks and Metropolis-Hastings Sampling in Modern Semantic Spaces

Human memory retrieval often resembles ecological foraging
where animals search for food in a “patchy” environment. Op-
timal foraging means strict adherences to the Marginal Value
Thereom (MVT) in which individuals exploit a “patch” of se-
mantically related concepts until it becomes less rewarding,
then switch to a new cluster. While human behavioral data
suggests foraging-like patterns in semantic fluency tasks, it
is still unknown whether modern high-dimensional embed-
ding spaces provide a sufficient representation for algorithms
to closely match observed human behavior. By leveraging
state-of-the-art embeddings and prior clustering and human se-
mantic fluency data I find that random walks on these seman-
tic embedding spaces produces results consistent with optimal
foraging and the MVT. Suprisingly, introducing Metroplois-
Hastings, an adaptive algorithm expected to model strategic
acceptance and rejection of new clusters, does not produce re-
sults consistent with observed human behavior. These findings
challenge the assumption that sophisticated sampling mecha-
nisms inherently provide better cognitive models of memory
retrieval. Instead, they highlight that appropriately structured
semantic embeddings, even with minimalistic sampling ap-
proaches, can produce near-optimal foraging dynamics. In do-
ing so, my results support the perspective of Hills (2012) rather
than Abbott (2015), demonstrating that modern embeddings
can approximate human memory foraging without relying on
complex acceptance criteria
