Clustering: Introduction
Clustering: Introduction
Clustering: Introduction
Author (@slack): Chioma Onyido (@Omabekee)
GitHub Repo: https://github.com/Omabekee/hackbio-cancer-internship
Clustering is a type of unsupervised learning which involves finding patterns within data distribution and regrouping the data into multiple clusters according to a defined distance measure so that the elements in one cluster are highly similar (intra-similarity) and the similarity between clusters (inter-similarity) is as tiny as possible (Oyewole & Thopil, 2023).


2392282286000










742883133116Figure 1: Clustering (Ezugwu et al., 2020)Figure 1: Clustering (Ezugwu et al., 2020)


Classification of clustering algorithms
There are five main categories of clustering algorithms: Partition-based clustering, Hierarchical clustering, Model-based clustering, Density-based clustering, and Grid-based clustering (Mehta et al., 2020).


41910010935400










right10962Figure 2: Classification of Clustering Algorithms (Mehta et al., 2020)Figure 2: Classification of Clustering Algorithms (Mehta et al., 2020)

Cluster analysis workflow
Figure 3 describes a typical workflow of clustering, beginning from data collection down to the validation of identified clusters using specific metrics such as the rand, silhouette, and Dunn metrics, to mention a few. Subsequently, bioinformatics analysis methods are applied to the identified clusters to annotate genes and pathways and identify potential relationships between those genes and specific biological events (Pineda et al., 2021). 
-49296000











502519285148Figure 3: Cluster Analysis Pipeline (Chioma Onyido)Figure 3: Cluster Analysis Pipeline (Chioma Onyido)


Applications of Clustering
Identification of tumour subtypes/ biomarkers for early diagnosis: Early diagnosis significantly improves patient survival, reduces treatment costs, and increases our understanding of the tumour immune microenvironment to inform more personalized immunotherapy approaches (Fahami et al., 2021). Various domains where clustering algorithms have been applied in this specific area include the detection of genes that have a significant impact on the mortality rate of colon cancer (Fahami et al., 2021) and identification of lung cancer subtypes using gene expression and clinical datasets (Kasa et al., 2020).
Identification of tumour subtypes for risk assessment and prediction: The heterogeneity of cancers among individuals is evidence that cancer subtypes exist. Clustering algorithms coupled with supervised machine learning algorithms are used to categorize patients into specific low-risk and high-risk groups and then propose prediction models or risk score calculators (Guo et al., 2021).
Image segmentation: K-means and hierarchical clustering have been employed to cluster functional images to identify high-risk sub-volumes in lung cancer ADDIN CSL_CITATION {"citationItems":[{"id":"ITEM-1","itemData":{"DOI":"https://doi.org/10.1016/j.radonc.2017.09.041","ISSN":"0167-8140","abstract":"Background and purpose\nWe aimed to identify tumour subregions with characteristic phenotypes based on pre-treatment multi-parametric functional imaging and correlate these subregions to treatment outcome. The subregions were created using imaging of metabolic activity (FDG-PET/CT), hypoxia (HX4-PET/CT) and tumour vasculature (DCE-CT).\nMaterials and methods\n36 non-small cell lung cancer (NSCLC) patients underwent functional imaging prior to radical radiotherapy. Kinetic analysis was performed on DCE-CT scans to acquire blood flow (BF) and volume (BV) maps. HX4-PET/CT and DCE-CT scans were non-rigidly co-registered to the planning FDG-PET/CT. Two clustering steps were performed on multi-parametric images: first to segment each tumour into homogeneous subregions (i.e. supervoxels) and second to group the supervoxels of all tumours into phenotypic clusters. Patients were split based on the absolute or relative volume of supervoxels in each cluster; overall survival was compared using a log-rank test.\nResults\nUnsupervised clustering of supervoxels yielded four independent clusters. One cluster (high hypoxia, high FDG, intermediate BF/BV) related to a high-risk tumour type: patients assigned to this cluster had significantly worse survival compared to patients not in this cluster (p = 0.035).\nConclusions\nWe designed a subregional analysis for multi-parametric imaging in NSCLC, and showed the potential of subregion classification as a biomarker for prognosis. This methodology allows for a comprehensive data-driven analysis of multi-parametric functional images.","author":[{"dropping-particle":"","family":"Even","given":"Aniek J G","non-dropping-particle":"","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Reymen","given":"Bart","non-dropping-particle":"","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Fontaine","given":"Matthew D","non-dropping-particle":"La","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Das","given":"Marco","non-dropping-particle":"","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Mottaghy","given":"Felix M","non-dropping-particle":"","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Belderbos","given":"José S A","non-dropping-particle":"","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Ruysscher","given":"Dirk","non-dropping-particle":"De","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Lambin","given":"Philippe","non-dropping-particle":"","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Elmpt","given":"Wouter","non-dropping-particle":"van","parse-names":false,"suffix":""}],"container-title":"Radiotherapy and Oncology","id":"ITEM-1","issue":"3","issued":{"date-parts":[["2017"]]},"page":"379-384","title":"Clustering of multi-parametric functional imaging to identify high-risk subvolumes in non-small cell lung cancer","type":"article-journal","volume":"125"},"uris":["http://www.mendeley.com/documents/?uuid=a69f5d4a-dcf4-434e-8410-3db8c0fc41f4"]}],"mendeley":{"formattedCitation":"(Even et al., 2017)","plainTextFormattedCitation":"(Even et al., 2017)","previouslyFormattedCitation":"(Even et al., 2017)"},"properties":{"noteIndex":0},"schema":"https://github.com/citation-style-language/schema/raw/master/csl-citation.json"}(Even et al., 2017). DBSCAN and deep learning have been used to develop a high-precision lung nodule detector to enable reidentification of previous nodules during lung cancer screening ADDIN CSL_CITATION {"citationItems":[{"id":"ITEM-1","itemData":{"DOI":"10.1186/s12880-021-00594-4","ISSN":"1471-2342 (Electronic)","PMID":"33836677","abstract":"BACKGROUND: Reidentification of prior nodules for temporal comparison is an  important but time-consuming step in lung cancer screening. We develop and evaluate an automated nodule detector that utilizes the axial-slice number of nodules found in radiology reports to generate high precision nodule predictions. METHODS: 888 CTs from Lung Nodule Analysis were used to train a 2-dimensional (2D) object detection neural network. A pipeline of 2D object detection, 3D unsupervised clustering, false positive reduction, and axial-slice numbers were used to generate nodule candidates. 47 CTs from the National Lung Cancer Screening Trial (NLST) were used for model evaluation. RESULTS: Our nodule detector achieved a precision of 0.962 at a recall of 0.573 on the NLST test set for any nodule. When adjusting for unintended nodule predictions, we achieved a precision of 0.931 at a recall 0.561, which corresponds to 0.06 false positives per CT. Error analysis revealed better detection of nodules with soft tissue attenuation compared to ground glass and undeterminable attenuation. Nodule margins, size, location, and patient demographics did not differ between correct and incorrect predictions. CONCLUSIONS: Utilization of axial-slice numbers from radiology reports allowed for development of a lung nodule detector with a low false positive rate compared to prior feature-engineering and machine learning approaches. This high precision nodule detector can reduce time spent on reidentification of prior nodules during lung cancer screening and can rapidly develop new institutional datasets to explore novel applications of computer vision in lung cancer imaging.","author":[{"dropping-particle":"","family":"Chillakuru","given":"Yeshwant Reddy","non-dropping-particle":"","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Kranen","given":"Kyle","non-dropping-particle":"","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Doppalapudi","given":"Vishnu","non-dropping-particle":"","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Xiong","given":"Zhangyuan","non-dropping-particle":"","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Fu","given":"Letian","non-dropping-particle":"","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Heydari","given":"Aarash","non-dropping-particle":"","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Sheth","given":"Aditya","non-dropping-particle":"","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Seo","given":"Youngho","non-dropping-particle":"","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Vu","given":"Thienkhai","non-dropping-particle":"","parse-names":false,"suffix":""},{"dropping-particle":"","family":"Sohn","given":"Jae Ho","non-dropping-particle":"","parse-names":false,"suffix":""}],"container-title":"BMC medical imaging","id":"ITEM-1","issue":"1","issued":{"date-parts":[["2021","4"]]},"language":"eng","page":"66","title":"High precision localization of pulmonary nodules on chest CT utilizing axial  slice number labels.","type":"article-journal","volume":"21"},"uris":["http://www.mendeley.com/documents/?uuid=4359b0bd-50cf-4957-9937-cbe69d46c9a9"]}],"mendeley":{"formattedCitation":"(Chillakuru et al., 2021)","plainTextFormattedCitation":"(Chillakuru et al., 2021)","previouslyFormattedCitation":"(Chillakuru et al., 2021)"},"properties":{"noteIndex":0},"schema":"https://github.com/citation-style-language/schema/raw/master/csl-citation.json"}(Chillakuru et al., 2021).

right6059200









right91774Figure 4: K-means Brain Tumour Example in R (Technophiles, 2020)Figure 4: K-means Brain Tumour Example in R (Technophiles, 2020)



Conclusion
Clustering algorithms are essential tools in cancer research. They help us pinpoint tumour subtypes, predict patient outcomes, classify the severity of disease complications using images and tailor treatments to individual needs. As we continue to explore these tools, they’re proving to be indispensable in pushing cancer research forward and improving patient care.

REFERENCES
Chillakuru, Y. R., Kranen, K., Doppalapudi, V., Xiong, Z., Fu, L., Heydari, A., Sheth, A., Seo, Y., Vu, T., & Sohn, J. H. (2021). High precision localization of pulmonary nodules on chest CT utilizing axial  slice number labels. BMC Medical Imaging, 21(1), 66. https://doi.org/10.1186/s12880-021-00594-4
Even, A. J. G., Reymen, B., La Fontaine, M. D., Das, M., Mottaghy, F. M., Belderbos, J. S. A., De Ruysscher, D., Lambin, P., & van Elmpt, W. (2017). Clustering of multi-parametric functional imaging to identify high-risk subvolumes in non-small cell lung cancer. Radiotherapy and Oncology, 125(3), 379–384. https://doi.org/https://doi.org/10.1016/j.radonc.2017.09.041
Ezugwu, A.E., Shukla, A.K., Agbaje, M.B., Oyelade, O.N., Jose-Garcia, A. $ Agushaka, J.O. (2021). Automatic clustering algorithms: a systematic review and bibliometric analysis of relevant literature. Neural Computing and Applications, 33, 6247-6306.
Fahami, M. A., Roshanzamir, M., Izadi, N. H., Keyvani, V., & Alizadehsani, R. (2021). Detection of effective genes in colon cancer: A machine learning approach. Informatics in Medicine Unlocked, 24, 100605. https://doi.org/https://doi.org/10.1016/j.imu.2021.100605
Guo, C., Wang, J., Wang, Y., Qu, X., Shi, Z., Meng, Y., Qiu, J., & Hua, K. (2021). Novel artificial intelligence machine learning approaches to precisely predict survival and site-specific recurrence in cervical cancer: A multi-institutional study. Translational Oncology, 14(5), 101032. https://doi.org/https://doi.org/10.1016/j.tranon.2021.101032
Kasa, S.R., Bhattacharaya, S. & Rajan, V. (2020). Gaussian mixture copulas for high-dimensional clustering and dependency-based subtyping. Bioinformatics, 36(2), 621-628.
Mehta, V., Bawa, S., & Singh, J. (2020). Analytical review of clustering techniques and proximity measures. Artificial Intelligence Review, 53, 5995 - 6023.
Oyewole, G.J. & Thopil, G.A. (2023). Data clustering: application and trends. Artificial Intelligence Review, 56, 6439–6475. https://doi.org/10.1007/s10462-022-10325-y 
Pineda, L.A., Pourshafeie, A., Ioannidis, A., Leibold, C. M., Chan, A. L., Bustamante, C. D., Frankovich, J., & Wojcik, G. L. (2021). Discovering prescription patterns in pediatric acute-onset neuropsychiatric syndrome patients. Journal of Biomedical Informatics, 113, 103664. https://doi.org/https://doi.org/10.1016/j.jbi.2020.103664
Technophiles (2020). K means clustering algorithm; K means brain tumor example in R and machine learning algorithms. Accessed from https://www.youtube.com/watch?reload=9&app=desktop&v=aqQrSNI8qGI on September 2nd, 2024.

 
