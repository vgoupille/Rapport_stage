
- Study DOL in PSR401 trough scRNA seq
- intership Valentin Goupille 
- Msc in bioinfo 

- ANR project divide 

-Division of Labour
	- notion econo, ecologique 
	- population s'orgarnise de maniere spontané ou non pour optimiser tache, spécialisation pour acomplir certainne tache , optimisation de fitness 
	- ex : travail,
	- optimisation 
	- peut s'observer a differente echelles 
		- echelle population : homme , fourmies (super-organisme : we can consider then à super antitie)
		- echelle organes chez eucaryotes 
		- echelles cellulaire : each cellules nerveruse, cellules musculaires , proccesus de differenciation permet adap, plantes aussi 
		- a echelle cellulaire stomates , ex chez les plantes C4  : compartimentation cellulaire 
	- microorganisme DOL (interspecies) => verif si rentre bien dans le concept
		-mettre figure division du travail
	- relation mutaliste sont elles DOL ? 
	- limites du concept restes floue



- mecanisme of DOL chez procaryotes restes encores peu etudier
	- encore moins au sein d'une seule espce
	- ANR projet DIvide collaboration between Stephane Hacquart and Philippes 
	- cherche a étudier ce phénomène chez pseudomonas braciss, voir si existe ou pas
- 
- 


figure du projet ANR 
- 3 WP : 
	- WP1 : nathan Vannier , they work 
	- WP2 : 
	- WP3 : they study DOL with 
	- WP2 : 

- Pourquoi cette souche ? 
	- colonisatrice racine d'arabido, drive le microbiote 
	- pathogène opportuniste 
	- contrairement majo  pas TSS3 
	- possède son genome 
	- possède une banque de mutants (rb )
- 3 genes majeurs  : qui drive cette fontion


- callout Attention dans ce WP : pas travaille en lien avec plantes , 
	- ici en milieu controlé 
		- voir si noisy regulation 
		- ou bien si certain cluster specialisé
- => pour faire ça neccesite d'un etudier les cellules  unique not bulk 
	- even if We work on isogenetique population 
	- acces Omics : dans le cadre mon stage we focus on transcrit  
		- Genomics, Transcriptics
		- in a first par we focus on transcripti =RNA 
		- we can imagine genetique mutation 
	- Bulk RNA-seq when you study all pop cell
	- Single cell when you want to study individual cell
		- => Single Cell RNA est vraiment beaucoup utiliser chez eucaryotes, (10X), Multi-omics, SC-ATAC seq , also begenin of sc proteomics metabolomics , spatial omics...
		- pour medecine notament et analyse devellopemental, and immunity 
	- IN procaryote vraiment loin de ça 
		- de ce que je sais scGenomic et transcriptomics
		- and it's not just a question of money but many technical challenge 
			- as Yann show you contamination for genomic
			- size cell, low RNA content , high rRNA content , cell wall, data analysis ...

		- 3  types methods exist to do 
			- each have is own advantage and inconvenient 
			- isolation and 
			- bac-seq 
			- split-pool => strat choice by solene
				- inconvenient big losing of cell
				- different step permeabilisation
				- 3 rounds 
					- barcode and UMI to identify 
					- combinaison => we can identify more than 
					-random hexamer and not random 
			- attention we amplifi all RNA not mRNA 
		- it's exist different methods depletion , but for une premiere approche pas choisi cette option, because we want result and 
		
			-contrary to bulk RNA-seq we don't have many read 

		- 1 milliard 5 de read visé , illumina sequencing 

- Differentes conditions :
	- fer / milieux restraint glucose 
		- why le fer : 
		- essential fonction for colonisation of racine ,
		- milieux restrreint autour des racines 
		- (milieux homogenes , agitation) 
	- suivi temporel (phase expo)
	- plusieurs replicat techni
	- plusieurs repli biologique 
- 
- voir si il y a noisy metabolis , ou pas 

each cell can't exist "alone" , just do limited function, need other cell to be "viable"
link to quorum sensing notion , => for that we need to have very hight depth sequencing 

we can also visualise very 






---

title: "Study of Division of Labor in Pseudomonas through single-cell RNA-seq"

subtitle: "From population-level to single-cell analysis"

author: "Valentin Goupille"

date: "2025-04-09"

format:

revealjs:

theme: simple

transition: slide

slide-number: true

footer: "Presentation of the internship"

css: custom.css

highlight-style: github

code-block-bg: true

code-block-border-left: "#31BAE9"

callout-appearance: simple

callout-icon: true

---



## General Context

- Division of Labour (DoL): concept from economy/sociology and biologie/ecology
     **Contexte scientifique** :
    
    La division du travail (DOL) est un concept central en économie, en écologie et en biologie. Il repose sur l’idée que des entités (individus, cellules, espèces…) se spécialisent dans des tâches spécifiques afin d’optimiser la performance collective ou la fitness globale.
    
    
- 
- ### Definition of Division of Labor (DOL)

- Specialization of tasks within a group

- Optimization of resource usage

- Enhanced collective performance

- Trade-offs between efficiency and autonomy

- => Optimization of task distribution to improve fitness and survival
    
- Observed at various levels:
    
    - Societal (humans, ants as super-organisms) (images de fourmis à droite)
        
    - **Organisme** : spécialisation des organes et des cellules  Tissue/organ level (e.g. eukaryotic organs) (images de tissus de brain et musculaire)
        
    - Cellular (e.g. C4 plants with bundle sheath cell and mesophyll cell) (images ) 
    - https://search-static.byjusweb.com/question-images/toppr_ext/questions/733605_677619_ans_e37b07c1068e4076a7489dc7a6760f48.jpg
        
    - Microbial communities (inter- and intra-species) (figure )
	- **Micro-organismes** :
        
        - DOL interspécifique (mutualisme, cross-feeding)
            
        - Reste à explorer la DOL **intraspécifique** chez les procaryotes.


:::



![[Article/image/giriDefiningDivisionLabor2019/image-8-x56-y237.png]]
[@giriDefiningDivisionLabor2019]  




::: notes


Le concept de Division du Travail trouve ses racines dans les travaux de plusieurs penseurs économiques importants :

1. Adam Smith (1776)

- Considéré comme le père fondateur de ce concept en économie

- Dans son ouvrage "La Richesse des Nations" (1776), il a développé l'idée que la division du travail augmente la productivité

- Il a illustré ce concept avec l'exemple célèbre de la manufacture d'épingles, montrant comment la spécialisation des tâches permet d'augmenter considérablement la production

Je vais vous expliquer pourquoi il y a un compromis (trade-off) entre l'efficacité et l'autonomie dans la Division du Travail.

Ce trade-off est un concept fondamental qui s'applique à plusieurs niveaux :

1. Au niveau économique

- Plus on se spécialise, plus on devient efficace dans une tâche spécifique

- Mais cette spécialisation réduit l'autonomie car on devient dépendant des autres pour les tâches qu'on ne maîtrise plus

- Exemple : Un artisan qui fabriquait tout son produit devient dépendant des fournisseurs de matières premières

### Eukaryotic Systems

  

- Well-documented phenomenon

- Classic example: cellular differentiation

- Specialized functions in multicellular organisms

- Social models (e.g., ants)

- Task optimization

- Super-organism organization

- Clear role distribution


1. Au niveau biologique (plus pertinent pour votre présentation)

- Dans les systèmes multicellulaires :

- Les cellules spécialisées sont très efficaces dans leur fonction spécifique

- Mais elles perdent leur autonomie et dépendent des autres cellules pour leur survie

- Exemple : Les cellules nerveuses sont très efficaces pour la transmission d'informations mais dépendent des cellules gliales pour leur survie

1. Dans les systèmes microbiens (votre sujet d'étude)

- Les bactéries spécialisées peuvent être très efficaces dans une fonction métabolique spécifique

- Mais cette spécialisation les rend dépendantes des autres membres de la communauté

- Exemple : Dans un biofilm, certaines bactéries peuvent se spécialiser dans la production de nutriments mais deviennent dépendantes des autres pour la protection ou la motilité

Ce trade-off est particulièrement important dans votre étude car :

- Il influence la stabilité des communautés microbiennes

- Il affecte la résilience face aux changements environnementaux

- Il détermine l'équilibre entre coopération et compétition dans la population

C'est pourquoi il est crucial de comprendre ce compromis dans votre analyse de la Division du Travail chez Pseudomonas brassicacearum.

diapoMM_Val.qmd

ants, termites 

Un **superorganisme** , ou **supraorganisme** , [[ 1 ]](https://en.wikipedia.org/wiki/Superorganism#cite_note-1) est un groupe d' organismes de la même [espèce qui interagissent](https://en.wikipedia.org/wiki/Species "Espèces") [de manière synergétique](https://en.wikipedia.org/wiki/Synergy "Synergie") . Une [communauté](https://en.wikipedia.org/wiki/Community_\(ecology\) "Communauté (écologie)") d'organismes d'espèces différentes qui interagissent de manière synergétique est appelée un _[holobionte](https://en.wikipedia.org/wiki/Holobiont "Holobiont")_ .

Le terme « superorganisme » est le plus souvent utilisé pour décrire une unité sociale d' animaux [eusociaux](https://en.wikipedia.org/wiki/Eusociality "Eusocialité") , où [la division du travail](https://en.wikipedia.org/wiki/Division_of_labour "Division du travail") est hautement spécialisée et où les individus ne peuvent survivre seuls pendant de longues périodes. [Les fourmis](https://en.wikipedia.org/wiki/Ants "Fourmis") sont l'exemple le plus connu d'un tel superorganisme. Un superorganisme peut être défini comme « un ensemble d'agents capables d'agir de concert pour produire des phénomènes régis par le collectif », [[ 2 ]](https://en.wikipedia.org/wiki/Superorganism#cite_note-2) phénomènes désignant toute activité « désirée par la ruche », comme la collecte de nourriture par les fourmis et [l'évitement des prédateurs](https://en.wikipedia.org/wiki/Antipredator_adaptation "Adaptation antiprédatrice") , [[ 3 ]](https://en.wikipedia.org/wiki/Superorganism#cite_note-3) [[ 4 ]](https://en.wikipedia.org/wiki/Superorganism#cite_note-4) ou le choix d'un nouveau site de nidification par les abeilles. [[ 5 ]](https://en.wikipedia.org/wiki/Superorganism#cite_note-5)

https://en.wikipedia.org/wiki/Holobiont

https://en.wikipedia.org/wiki/Superorganism#:~:text=7%20External%20links-,Concept,example%20of%20such%20a%20superorganism.


https://pmc.ncbi.nlm.nih.gov/articles/PMC2948986/#:~:text=The%20division%20of%20labor%20between,complexes%20and%20other%20metabolic%20pathways.


Un **holobionte** est un assemblage d'un [hôte](https://en.wikipedia.org/wiki/Host_\(biology\) "Hôte (biologie)") et des nombreuses autres espèces vivant à l'intérieur ou autour de lui, qui forment ensemble une [unité écologique](https://en.wikipedia.org/wiki/Ecological_unit "Unité écologique") discrète par [symbiose](https://en.wikipedia.org/wiki/Symbiosis "Symbiose") , [[ 2 ]](https://en.wikipedia.org/wiki/Holobiont#cite_note-Margulis1991-2) bien qu'il existe une controverse sur cette discrétion. Les composants d'un holobionte sont des espèces individuelles ou [bionts , tandis que le](https://en.wiktionary.org/wiki/biont "wikt:biont") [génome](https://en.wikipedia.org/wiki/Genome "génome") combiné de tous les bionts est l' [hologénome](https://en.wikipedia.org/wiki/Hologenome_theory_of_evolution "Théorie de l'évolution par l'hologénome") .

link to symbiotic relation: mutualist, parasitism, 
https://en.wikipedia.org/wiki/Biological_interaction

aussi mitochondries, chloroplastes 
::: 


Roots of healthy plants host diverse bacteria that are collectively referred to as the bacterial root microbiota. Unlike their multicellular eukaryotic hosts that evolved diverse cell-types to achieve distinct biological functions and promote a division of labour, unicellular organisms such as bacteria rely on metabolic exchange(s) with their surrounding biotic environment. 

Recent reports, including our own work (Mataigne et al. 2021, Mataigne et al. 2022), indicate that metabolic interdependencies and cross feeding exchanges are widespread among taxonomically diverse bacteria and likely drive microbial co-existence within complex bacterial communities (e.g. Estrella et al., 2016; Adkins-Jablonsky et al., 2021).


![[Pasted image 20250408091911.png]]



“However, a major unsolved question is whether populations of genetically identical bacteria can minimise energetically costly processes by each executing different metabolic tasks at the intra-population level. Here, we hypothesise that metabolic cooperation between bacterial intra-populations plays a key role in modulating population dynamics, competitiveness and persistence at the root soil-interface. This hypothesis also builds on the idea that within a population, bacteria are inclined to 'noisy regulation' of metabolism (Lopez & Wingreen, 2022), i.e. they do not all adjust their genome expression to the environmental constraints in the same way. In a stable environment, this would be expected to limit bacterial growth. However, because bacteria excrete compounds (‘leaky function’ forming a ‘metabolic marketplace’), the selection pressures at the population level would favour intra-population cross-feeding(s) fitting to the constraints, either supported by variations in gene expression or/and selection for useful variants that arise from non-synonymous mutations (Figure 1). Thus, by ‘noise-averaging cooperation’ (Lopez & Wingreen, 2022) an intra-population-level division of labour is expected.” 


## Concept Overview

  

### Traditional View

  

- All bacterial cells identical

- Clonal populations

- Uniform behavior

  

### New Paradigm

  

- Population heterogeneity

- Metabolic specialization

- Task distribution

- Community-level benefits
une diapo sur les differents mecanismes : 

### **Projet ANR DIVIDE**

- **Objectif global** : Étudier la DOL au sein d’une population isogénique de _Pseudomonas brassicacearum_ (PsR401)
    
- **Partenariat** : Stéphane Hacquard (MPIPZ) 
    
- **Trois workpackages (WP)** : with Solène and Philippe we focus on WP2 
    
differentes approches focus sur niveaux transcriptomiques


![[Pasted image 20250408094156.png]]




# 🌱 Microbial Strategies for Population-Level Resilience

| Concept                | Description                                                                 | Example Behavior                         | Key Feature                            |
|------------------------|-----------------------------------------------------------------------------|------------------------------------------|-----------------------------------------|
| 🧩 **Metabolic Cooperation** | Cells **specialize** in different metabolic tasks to benefit the group.       | One cell produces siderophores, another breaks down sugars. | **Functional division of labor**       |
| 🎲 **Bet-Hedging**           | Random expression of stress-related genes allows **some cells** to survive unpredictable changes. | A few cells express stress genes before stress appears.      | **Phenotypic diversification**         |
| 🔀 **Noisy Regulation**      | Gene expression **varies stochastically** between identical cells.             | Random fluctuations in gene activity, not coordinated.       | **Uncoordinated heterogeneity**        |


## 🧬 Population of Genetically Identical Cells


| Cell A      | Cell B      | Cell C  |
| ----------- | ----------- | ------- |
| Iron Uptake | Sugar Util. | Dormant |
| Genes ↑     | Genes ↑     |         |


➡️ **Metabolic cooperation**: specialization in metabolic functions (e.g., iron vs. sugar).  
➡️ **Noisy regulation**: spontaneous entry into a dormant or stress-ready state.  
➡️ All cells are genetically identical → **phenotypic diversity** increases survival.





However, a major unsolved question is whether populations of genetically identical bacteria can minimise energetically costly processes by each executing different metabolic tasks at the intra-population level. Here, we hypothesise that metabolic cooperation between bacterial intra-populations plays a key role in modulating population dynamics, competitiveness and persistence at the root soil-interface. This hypothesis also builds on the idea that within a population, bacteria are inclined to 'noisy regulation' of metabolism (Lopez & Wingreen, 2022), i.e. they do not all adjust their genome expression to the environmental constraints in the same way. In a stable environment, this would be expected to limit bacterial growth. However, because bacteria excrete compounds (‘leaky function’ forming a ‘metabolic marketplace’), the selection pressures at the population level would favour intra-population cross-feeding(s) fitting to the constraints, either supported by variations in gene expression or/and selection for useful variants that arise from non-synonymous mutations (Figure 1). Thus, by ‘noise-averaging cooperation’ (Lopez & Wingreen, 2022) an intra-population-level division of labour is expected.


limites floue entre ces differents concepts : et autres concepts. : spatio-temporel, sinetique de changement , epigentique , phenotypique switch ..., quorum sensing 
(sRNA, integrase,  ) https://fr.wikipedia.org/wiki/Intégron 
https://fr.wikipedia.org/wiki/Recombinase

autre article :  [@dekkersSitespecificRecombinaseRequired1998]
lien avec [[Pseudomonas]] , [[integrase]]

[[article intergrase pseudomonas variation phase ]]
[[integrase]]


-etudier en utilisant scRNA-seq 

- Funded research program investigating Division of Labour in prokaryotes
    
- Model organism: _Pseudomonas brassicacearum_ (PSR401)
    
    - Root colonizer of _Arabidopsis_
        
    - Lacks T3SS, unlike many _Pseudomonas_
        
    - Genetically tractable: available RB mutant library (WP1)



### **Cas d’étude : Pseudomonas brassicacearum R401 (PsR401)**

  

Cette souche bactérienne, connue pour coloniser efficacement les racines d’_Arabidopsis thaliana_, utilise **plusieurs stratégies indépendantes** et **complémentaires** pour :

1. **Entrer en compétition avec d’autres microbes**, via :
    
    - Production d’un **antimicrobien**.
        
    - Production d’une **molécule qui capte le fer**, un nutriment essentiel mais limité dans la rhizosphère.
        
    
2. **Coloniser la racine et provoquer une certaine pathogénicité**, via :
    
    - Production d’un **phytotoxine** (non encore publiée) qui agit à la fois sur la plante hôte et sur l’efficacité de la colonisation.
        
    

  

💡 **Ces mécanismes fonctionnent de manière indépendante mais additive** : chacun contribue à la colonisation, et leur combinaison renforce la persistance de la souche à la racine.



### Short description of Pseudomonas brassicacearum R401

faire recherche sur pseudomonas brassicacearum R401
on possède son genome (ici il s'agit de son genome issu d'un assemblage )
=> discussion avec Solène : évolution rapide peut etre ne voit pas toutes les spécificité 


article de Felix : sur pseudomonas braci selection des genes 
[@getzkeCofunctioningBacterialExometabolites2023]
#Important 

R401 est très concurrentiel (grace à ses 3 facteurs : toxine, chelateur, ...) (felix)

R401 est isolé de plantes saine (le microbiote avec autre pseudomonas protege la plante)

R401 est opportunistes (en condition de sel devient pathogène grace notamment à sa toxine brassicapeptin : exomethabolite inhibiteur qui troue les parois (vegetale et aussi autres bacteries ) ) 


- [@mesnyCoevolutionPlantHolobiont2023]parle rapidement de pseudomonas ::
	
- “A similar analysis of binary interactions among 198 A. thaliana root-isolated bacteria identified Pseudomonas brassicacearum as a potent antagonist that relies on the combined action of two exometabolites, an antimicrobial and an iron chelator, that suppress competitors, and thereby promote its root colonization (Getzke et al, 2023).” (Mesny et al., 2023, p. 9)
- “Resource competition represents an important mechanism for indirect microbial antagonism (Hassani et al, 2018). The ability to rapidly utilize a limited resource can be detrimental for a less competitive microbe. Additionally, microbes can sequester resources, preventing utilization by other community members. Siderophores chelate soil iron for microbial uptake and are differentially exploitable to prevent other microbes to obtain essential iron resources (Joshi et al, 2006).” (Mesny et al., 2023, p. 9)
- “For instance, pyoverdine is an iron-chelator produced by some A. thaliana root-associated bacteria contributing to their competitiveness (Getzke et al, 2023).” (Mesny et al., 2023, p. 9)










- la souche pseudomonas R401 peut etre pathogène en condition de sel ; osmolyte ... : [@getzkePhysiochemicalInteractionOsmotic2024]
-  Unlike many pathogenic bacteria, R401 lacks genes for a type III secretion system and does not overgrow or suppress plant immune responses.
- 
• Research by Stéphane Hacquard (Max Planck Institute for Plant Breeding Research, Cologne) found that ~95% of plant-associated bacteria are neutral or beneficial in one-on-one interactions with _Arabidopsis thaliana_.

• _Pseudomonas brassicacearum_ R401 is a dominant plant microbiota member but can be detrimental under laboratory conditions.

• However, in natural soil conditions, R401 does not cause disease, suggesting it requires specific environmental conditions to become pathogenic.

• This mutation turned R401 into a plant-beneficial bacterium.

• **Effects of brassicapeptin**

• Alone, brassicapeptin can induce plant disease under high salt conditions.

• It is toxic to _Arabidopsis_, tomato plants under salt stress, and other microbes.

• The molecule has a fatty acid tail linked to amino acids and can form pores in plant membranes, likely increasing salt stress sensitivity.



Successful establishment of bacteria at roots is driven by multiple independent biological processes involved in both host-microbe (i.e., signal recognition, chemotaxis, surface attachment, biofilm formation, virulence factors) and microbe-microbe interactions (i.e., production of antimicrobials or public goods) (Knights et al. 2021). Therefore, we postulate that the simultaneous activation of these diverse processes is costly and that cooperation between genetically identical strains is key for promoting bacterial pervasiveness at roots. Consistent with this hypothesis, we recently demonstrated that the robust root coloniser Pseudomonas brassicacearum R401 (hereafter referred to as PsR401) deploys multiple independent strategies that co-function to promote colonisation and persistence at roots (Getzke et al., 2023). Notably, we identified two independent processes involved in microbe-microbe competition – namely the production of an antimicrobial and of a molecule scavenging the micronutrient iron – that act additively to promote strain competitiveness at roots (Getzke et al., 2023). We identified a third bacterial locus in PsR401 involved in the production of a phytotoxin that promotes both pathogenicity and root colonisation in mono-association experiments with A. thaliana (unpublished, see below). Our work provides proof-of-concept data illustrating that PsR401 deploys diverse exo-metabolites during root microbiota establishment that have versatile functions in host-microbe and microbe-microbe interactions. Together with earlier work (Gu et al. 2020, Harbort et al 2020), it also delineates iron as a major micronutrient modulating strain competitiveness and proliferation at roots. Given that the public good iron becomes rate-limiting in the root compartment and that production of the above-mentioned processes are all modulated by iron availability (Lim et al., 2012; Palma et al., 2003, Mo et al. 1991), we anticipate that division of labour between bacterial intra-populations is bolstered under iron limiting conditions such as those found in the root habitat.

---

### **🧪** 

### **Fer : un facteur limitant clé dans la rhizosphère**

  

Le fer est un **“bien public”** : toutes les bactéries en ont besoin, mais il est rare dans l’environnement racinaire. Or :

- Plusieurs des mécanismes cités (production d’antimicrobiens, de phytotoxines, etc.) sont **régulés par la disponibilité en fer**.
    
- Dans des conditions de **carence en fer** (ce qui est fréquent au niveau des racines), **produire toutes ces fonctions à la fois devient encore plus coûteux**.
    

---

### **🤝** 

### **Lien avec la division du travail bactérienne**

  

Dans ce contexte, **une division du travail au sein de la population de PsR401** serait **avantageuse** :

- Certaines cellules se spécialiseraient dans la **production d’antimicrobiens**.
    
- D’autres dans la **capture du fer**.
    
- D’autres encore dans la **production de phytotoxines**.
    

  

➡️ Cela permettrait à la **population globale d’assurer toutes les fonctions nécessaires** à une colonisation efficace, **tout en minimisant le coût pour chaque cellule individuellement**.



# Scientific Questions

H2 (WP2): In a population of Pseudomonas brassicacearum (PsR401), individuals activate different biological and metabolic processes, promoting intra-population functional diversity and allowing a rapid division of labour to colonise the environment (with a particular focus on the root-environment)

- Is there functional specialization (DoL) within clonal populations of PSR401?
    
- Can scRNA-seq reveal transcriptional heterogeneity?
    
- Are certain cells specialized under stress (e.g. iron or glucose limitation)?
    
- Is this driven by noisy gene expression or regulated mechanisms?

The major goal of the project is to obtain mechanistic insights into whether and how a robust bacterial root coloniser can partition tasks within intra-populations to achieve greater functional diversity, minimise the energetic cost of processes, and adapt to environmental constraints to promote population density at roots. Given that at a micro-scale, niche partitioning is often observed between bacterial populations at roots, we hypothesise that cooperation within a bacterial population is at least as important as cooperation between bacterial populations for promoting overall population density. We anticipate that a single cell can only achieve a limited number of functions and that greater functional diversity can be achieved through within-population cooperation, especially under nutritional constraints. Given that iron is rate-limiting in roots, we particularly aim at testing the relevance of division of labour between bacterial intra-populations in the context of iron availability.





# Why Single-cell RNA-seq?

- Traditional bulk RNA-seq masks cell-to-cell variability
    
- scRNA-seq enables:
    
    - Identification of cell states and subpopulations
        
    - Characterization of heterogeneous responses
        

1. scRNA-seq in bacteria: methods, limitations
    

## Challenges in Bacteria

- No polyA tail on mRNA
    
- Very low RNA content per cell
    
- High rRNA content
    
- Cell wall makes lysis difficult
    

- 3 main methods: Isolation-based, BacDrop, Split-seq
    
- Focus on Split-seq and MicroSplit (Kuchina et al.)


# Experimental Setup

## Conditions Studied

- Iron limitation and glucose limitation
    
- Controlled homogeneous environment (no plant interaction)
    
- Time course during exponential phase
    
- Multiple biological and technical replicates
    

## Sequencing

- Split-seq protocol adapted to PSR401 (based on Kuchina et al.)
    
- Sequencing depth: 1.5 billion reads (Illumina NovaSeq)




		- in a first par we focus on transcripti =RNA 
		- we can imagine genetique mutation 
	- Bulk RNA-seq when you study all pop cell
	- Single cell when you want to study individual cell
		- => Single Cell RNA est vraiment beaucoup utiliser chez eucaryotes, (10X), Multi-omics, SC-ATAC seq , also begenin of sc proteomics metabolomics , spatial omics...
