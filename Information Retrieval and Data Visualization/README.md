**1. Information Retrieval:**

**Project: Hybrid Recommender System for Movie Recommendation with Explicit Negative Feedback**

This project focuses on the design and implementation of a hybrid recommender system for movies that combines Collaborative Filtering (CF) and Content-Based Filtering (CB) to provide robust and personalized recommendations. The collaborative component is based on an Alternating Least Squares (ALS) matrix factorization model trained on a sparse user–item interaction matrix, allowing the system to learn latent representations of users and movies from implicit feedback.

The content-based component leverages rich movie metadata extracted from a relational database, including textual descriptions, titles, directors, cast members, genres, release year, and duration. Textual features are encoded using TF-IDF vectorization and compared via cosine similarity, while categorical attributes are evaluated using Jaccard similarity. Numerical attributes such as release year and duration are incorporated through heuristic similarity functions. These heterogeneous similarities are combined into a single weighted similarity score, enabling fine-grained control over the contribution of each feature.

A key aspect of the system is the explicit handling of negative feedback: user preferences are modeled in binary form, and movies similar to those disliked by a user are penalized in the recommendation score. The final recommendation score is obtained through an adaptive hybrid strategy that blends normalized CF and CB scores, while automatically falling back to pure content-based recommendations in cold-start scenarios where collaborative information is insufficient.

The system is evaluated using synthetic users with controlled preference patterns, such as users favoring a specific movie franchise, to assess whether the recommender correctly infers missing preferences and suggests coherent related content. An interactive command-line interface is provided to explore recommendations, inspect scoring components, and analyze execution performance.

**Software & Tools:** Python, Pandas, NumPy, SciPy, scikit-learn, implicit, SQLite

**Models & Techniques Used:** Hybrid recommender systems; Collaborative Filtering with ALS; implicit feedback modeling; sparse user–item matrices; TF-IDF vectorization; cosine similarity; Jaccard similarity; feature-weighted similarity aggregation; cold-start detection; negative preference penalization; interactive recommendation evaluation.

**2. Data Visualization:**

**Project: Analysis of Drug Use and Fentanyl Enforcement**
This project focuses on the construction of a full data‑visualization pipeline to explore illegal drug‑use trends and fentanyl‑trafficking patterns in the United States. The work begins with the preprocessing and integration of two heterogeneous datasets: the 2021–2022 National Surveys on Drug Use and Health (NSDUH) and the HIDTA fentanyl‑seizure records (2017–2023). The data are cleaned, normalized, and structured to enable comparative analysis across age groups, substances, states, and regions.

The visualization component includes the design of multiple chart types tailored to different analytical questions. Bar charts and grouped comparisons are used to identify the most consumed substances by age group. Choropleth maps display geographic variation in drug use across U.S. states, using color hue and shade to encode prevalence. Small‑multiple maps allow cross‑substance comparison while preserving geographic consistency. Time‑series area plots illustrate the evolution of fentanyl seizures by region and by form (powder vs. pills), highlighting the rapid escalation in the Western United States.

The project emphasizes the selection and justification of visual channels—position, color, area, containment—to ensure clarity, comparability, and interpretability. The resulting visualizations support the identification of demographic patterns, regional disparities, and emerging trafficking trends.

**Software & Tools:** Python, Pandas, NumPy, Matplotlib, Seaborn, GeoPandas, Jupyter Notebook; JSON/GeoJSON for geographic boundaries.

**Models & Techniques Used:**
Data preprocessing and cleaning; exploratory data analysis; geospatial visualization (choropleth maps); time‑series visualization; small‑multiple design; categorical and quantitative encoding; color‑hue and color‑shade mapping; ranking and comparative analysis; visual‑channel optimization for effective communication.
