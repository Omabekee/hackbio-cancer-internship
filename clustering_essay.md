<!--StartFragment-->

# Decoding Cancer Patterns with Clustering

### Author (@slack): Chioma Onyido (@Omabekee)

Clustering is a type of unsupervised learning which involves finding patterns within data distribution and regrouping the data into multiple clusters according to a defined distance measure so that the elements in one cluster are highly similar (intra-similarity) and the similarity between clusters (inter-similarity) is as tiny as possible (Oyewole & Thopil, 2023).

![Figure 1: Clustering (Ezugwu _et al_., 2020)](images/clusterinng.png)

## Classification of clustering algorithms

There are five main categories of clustering algorithms: Partition-based clustering, Hierarchical clustering, Model-based clustering, Density-based clustering, and Grid-based clustering (Mehta _et al_., 2020).

![Figure 2: Classification of Clustering Algorithms (Mehta _et al_., 2020)](images/classification_of_clustering.png)


## Cluster analysis workflow

Figure 3 describes a typical workflow of clustering, beginning from data collection down to the validation of identified clusters using specific metrics such as the rand, silhouette, and Dunn metrics, to mention a few. Subsequently, bioinformatics analysis methods are applied to the identified clusters to annotate genes and pathways and identify potential relationships between those genes and specific biological events (Pineda _et al_., 2021). 

![Figure 3: Cluster Analysis Pipeline (Chioma Onyido)](images/cluster_analysis_pipeline.png) 

## Applications of Clustering

1. **Identification of tumour subtypes/ biomarkers for early diagnosis:** Early diagnosis significantly improves patient survival, reduces treatment costs, and increases our understanding of the tumour immune microenvironment to inform more personalized immunotherapy approaches (Fahami _et al_., 2021)**.** Various domains where clustering algorithms have been applied in this specific area include the detection of genes that have a significant impact on the mortality rate of colon cancer (Fahami _et al_., 2021) and identification of lung cancer subtypes using gene expression and clinical datasets (Kasa _et al_., 2020).

2. **Identification of tumour subtypes for risk assessment and prediction:** The heterogeneity of cancers among individuals is evidence that cancer subtypes exist. Clustering algorithms coupled with supervised machine learning algorithms are used to categorize patients into specific low-risk and high-risk groups and then propose prediction models or risk score calculators (Guo _et al_., 2021).

3. **Image segmentation:** K-means and hierarchical clustering have been employed to cluster functional images to identify high-risk sub-volumes in lung cancer (Even _et al_., 2017). DBSCAN and deep learning have been used to develop a high-precision lung nodule detector to enable reidentification of previous nodules during lung cancer screening (Chillakuru _et al_., 2021).

![Figure 4: K-means Brain Tumour Example in R (Technophiles, 2020)](images/kmeans_brain_tumor.jpg)

**Conclusion**

Clustering algorithms are essential tools in cancer research. They help us pinpoint tumour subtypes, predict patient outcomes, classify the severity of disease complications using images and tailor treatments to individual needs. As we continue to explore these tools, they’re proving to be indispensable in pushing cancer research forward and improving patient care.

**REFERENCES**

Chillakuru, Y. R., Kranen, K., Doppalapudi, V., Xiong, Z., Fu, L., Heydari, A., Sheth, A., Seo, Y., Vu, T., & Sohn, J. H. (2021). High precision localization of pulmonary nodules on chest CT utilizing axial  slice number labels. _BMC Medical Imaging_, _21_(1), 66. https\://doi.org/10.1186/s12880-021-00594-4

Even, A. J. G., Reymen, B., La Fontaine, M. D., Das, M., Mottaghy, F. M., Belderbos, J. S. A., De Ruysscher, D., Lambin, P., & van Elmpt, W. (2017). Clustering of multi-parametric functional imaging to identify high-risk subvolumes in non-small cell lung cancer. _Radiotherapy and Oncology_, _125_(3), 379–384. https\://doi.org/https\://doi.org/10.1016/j.radonc.2017.09.041

Ezugwu, A.E., Shukla, A.K., Agbaje, M.B., Oyelade, O.N., Jose-Garcia, A. $ Agushaka, J.O. (2021). Automatic clustering algorithms: a systematic review and bibliometric analysis of relevant literature. _Neural Computing and Applications_, _33_, 6247-6306.

Fahami, M. A., Roshanzamir, M., Izadi, N. H., Keyvani, V., & Alizadehsani, R. (2021). Detection of effective genes in colon cancer: A machine learning approach. _Informatics in Medicine Unlocked_, _24_, 100605. https\://doi.org/https\://doi.org/10.1016/j.imu.2021.100605

Guo, C., Wang, J., Wang, Y., Qu, X., Shi, Z., Meng, Y., Qiu, J., & Hua, K. (2021). Novel artificial intelligence machine learning approaches to precisely predict survival and site-specific recurrence in cervical cancer: A multi-institutional study. _Translational Oncology_, _14_(5), 101032. https\://doi.org/https\://doi.org/10.1016/j.tranon.2021.101032

Kasa, S.R., Bhattacharaya, S. & Rajan, V. (2020). Gaussian mixture copulas for high-dimensional clustering and dependency-based subtyping. _Bioinformatics_, _36_(2), 621-628.

Mehta, V., Bawa, S., & Singh, J. (2020). Analytical review of clustering techniques and proximity measures. _Artificial Intelligence Review, 53_, 5995 - 6023.

Oyewole, G.J. & Thopil, G.A. (2023). Data clustering: application and trends. _Artificial Intelligence Review, 56_, 6439–6475. https\://doi.org/10.1007/s10462-022-10325-y 

Pineda, L.A., Pourshafeie, A., Ioannidis, A., Leibold, C. M., Chan, A. L., Bustamante, C. D., Frankovich, J., & Wojcik, G. L. (2021). Discovering prescription patterns in pediatric acute-onset neuropsychiatric syndrome patients. _Journal of Biomedical Informatics_, _113_, 103664. https\://doi.org/https\://doi.org/10.1016/j.jbi.2020.103664

Technophiles (2020). K means clustering algorithm; K means brain tumor example in R and machine learning algorithms. Accessed from <https://www.youtube.com/watch?reload=9&app=desktop&v=aqQrSNI8qGI> on September 2<sup>nd</sup>, 2024.

 

\


<!--EndFragment-->
