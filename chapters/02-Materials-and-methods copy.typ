// Some definitions presupposed by pandoc's typst output.
#let blockquote(body) = [
  #set text( size: 0.92em )
  #block(inset: (left: 1.5em, top: 0.2em, bottom: 0.2em))[#body]
]

#let horizontalrule = [
  #line(start: (25%,0%), end: (75%,0%))
]

#let endnote(num, contents) = [
  #stack(dir: ltr, spacing: 3pt, super[#num], contents)
]

#show terms: it => {
  it.children
    .map(child => [
      #strong[#child.term]
      #block(inset: (left: 1.5em, top: -0.4em))[#child.description]
      ])
    .join()
}

// Some quarto-specific definitions.

#show raw.where(block: true): set block(
    fill: luma(230),
    width: 100%,
    inset: 8pt,
    radius: 2pt
  )

#let block_with_new_content(old_block, new_content) = {
  let d = (:)
  let fields = old_block.fields()
  fields.remove("body")
  if fields.at("below", default: none) != none {
    // TODO: this is a hack because below is a "synthesized element"
    // according to the experts in the typst discord...
    fields.below = fields.below.amount
  }
  return block.with(..fields)(new_content)
}

#let empty(v) = {
  if type(v) == "string" {
    // two dollar signs here because we're technically inside
    // a Pandoc template :grimace:
    v.matches(regex("^\\s*$")).at(0, default: none) != none
  } else if type(v) == "content" {
    if v.at("text", default: none) != none {
      return empty(v.text)
    }
    for child in v.at("children", default: ()) {
      if not empty(child) {
        return false
      }
    }
    return true
  }

}

// Subfloats
// This is a technique that we adapted from https://github.com/tingerrr/subpar/
#let quartosubfloatcounter = counter("quartosubfloatcounter")

#let quarto_super(
  kind: str,
  caption: none,
  label: none,
  supplement: str,
  position: none,
  subrefnumbering: "1a",
  subcapnumbering: "(a)",
  body,
) = {
  context {
    let figcounter = counter(figure.where(kind: kind))
    let n-super = figcounter.get().first() + 1
    set figure.caption(position: position)
    [#figure(
      kind: kind,
      supplement: supplement,
      caption: caption,
      {
        show figure.where(kind: kind): set figure(numbering: _ => numbering(subrefnumbering, n-super, quartosubfloatcounter.get().first() + 1))
        show figure.where(kind: kind): set figure.caption(position: position)

        show figure: it => {
          let num = numbering(subcapnumbering, n-super, quartosubfloatcounter.get().first() + 1)
          show figure.caption: it => {
            num.slice(2) // I don't understand why the numbering contains output that it really shouldn't, but this fixes it shrug?
            [ ]
            it.body
          }

          quartosubfloatcounter.step()
          it
          counter(figure.where(kind: it.kind)).update(n => n - 1)
        }

        quartosubfloatcounter.update(0)
        body
      }
    )#label]
  }
}

// callout rendering
// this is a figure show rule because callouts are crossreferenceable
#show figure: it => {
  if type(it.kind) != "string" {
    return it
  }
  let kind_match = it.kind.matches(regex("^quarto-callout-(.*)")).at(0, default: none)
  if kind_match == none {
    return it
  }
  let kind = kind_match.captures.at(0, default: "other")
  kind = upper(kind.first()) + kind.slice(1)
  // now we pull apart the callout and reassemble it with the crossref name and counter

  // when we cleanup pandoc's emitted code to avoid spaces this will have to change
  let old_callout = it.body.children.at(1).body.children.at(1)
  let old_title_block = old_callout.body.children.at(0)
  let old_title = old_title_block.body.body.children.at(2)

  // TODO use custom separator if available
  let new_title = if empty(old_title) {
    [#kind #it.counter.display()]
  } else {
    [#kind #it.counter.display(): #old_title]
  }

  let new_title_block = block_with_new_content(
    old_title_block, 
    block_with_new_content(
      old_title_block.body, 
      old_title_block.body.body.children.at(0) +
      old_title_block.body.body.children.at(1) +
      new_title))

  block_with_new_content(old_callout,
    block(below: 0pt, new_title_block) +
    old_callout.body.children.at(1))
}

// 2023-10-09: #fa-icon("fa-info") is not working, so we'll eval "#fa-info()" instead
#let callout(body: [], title: "Callout", background_color: rgb("#dddddd"), icon: none, icon_color: black) = {
  block(
    breakable: false, 
    fill: background_color, 
    stroke: (paint: icon_color, thickness: 0.5pt, cap: "round"), 
    width: 100%, 
    radius: 2pt,
    block(
      inset: 1pt,
      width: 100%, 
      below: 0pt, 
      block(
        fill: background_color, 
        width: 100%, 
        inset: 8pt)[#text(icon_color, weight: 900)[#icon] #title]) +
      if(body != []){
        block(
          inset: 1pt, 
          width: 100%, 
          block(fill: white, width: 100%, inset: 8pt, body))
      }
    )
}



#let article(
  title: none,
  subtitle: none,
  authors: none,
  date: none,
  abstract: none,
  abstract-title: none,
  cols: 1,
  margin: (x: 1.25in, y: 1.25in),
  paper: "us-letter",
  lang: "en",
  region: "US",
  font: "linux libertine",
  fontsize: 11pt,
  title-size: 1.5em,
  subtitle-size: 1.25em,
  heading-family: "linux libertine",
  heading-weight: "bold",
  heading-style: "normal",
  heading-color: black,
  heading-line-height: 0.65em,
  sectionnumbering: none,
  toc: false,
  toc_title: none,
  toc_depth: none,
  toc_indent: 1.5em,
  doc,
) = {
  set page(
    paper: paper,
    margin: margin,
    numbering: "1",
  )
  set par(justify: true)
  set text(lang: lang,
           region: region,
           font: font,
           size: fontsize)
  set heading(numbering: sectionnumbering)
  if title != none {
    align(center)[#block(inset: 2em)[
      #set par(leading: heading-line-height)
      #if (heading-family != none or heading-weight != "bold" or heading-style != "normal"
           or heading-color != black or heading-decoration == "underline"
           or heading-background-color != none) {
        set text(font: heading-family, weight: heading-weight, style: heading-style, fill: heading-color)
        text(size: title-size)[#title]
        if subtitle != none {
          parbreak()
          text(size: subtitle-size)[#subtitle]
        }
      } else {
        text(weight: "bold", size: title-size)[#title]
        if subtitle != none {
          parbreak()
          text(weight: "bold", size: subtitle-size)[#subtitle]
        }
      }
    ]]
  }

  if authors != none {
    let count = authors.len()
    let ncols = calc.min(count, 3)
    grid(
      columns: (1fr,) * ncols,
      row-gutter: 1.5em,
      ..authors.map(author =>
          align(center)[
            #author.name \
            #author.affiliation \
            #author.email
          ]
      )
    )
  }

  if date != none {
    align(center)[#block(inset: 1em)[
      #date
    ]]
  }

  if abstract != none {
    block(inset: 2em)[
    #text(weight: "semibold")[#abstract-title] #h(1em) #abstract
    ]
  }

  if toc {
    let title = if toc_title == none {
      auto
    } else {
      toc_title
    }
    block(above: 0em, below: 2em)[
    #outline(
      title: toc_title,
      depth: toc_depth,
      indent: toc_indent
    );
    ]
  }

  if cols == 1 {
    doc
  } else {
    columns(cols, doc)
  }
}

#set table(
  inset: 6pt,
  stroke: none
)

#show: doc => article(
  toc_title: [Table of contents],
  toc_depth: 3,
  cols: 1,
  doc,
)

= Materials and Methods
<sec-materials-and-methods>
z

== Bacterial culture
<bacterial-culture>
\-boite petri population isogenique , puis culture en liquid medium

#emph[Pseudomonas brassicacearum] R401 was grown in liquid medium…. at …°C.

Two different conditions were applied to the bacteria :

- Low glucose and low iron (M9 medium)

- High glucose and high iron (M9F medium)

For each condition, 3 replicates were grown (biological replicates) and des cellules ont été prelevées de ces cultures. The DO was measured during the growth à 3 DO\_timepoints (OD 0.1, 0.2, 0.3) Which do a total of 18 biological Samples/ Conditions (2 medium \* 3 biological replicates \* 3 timepoints)

Culture medium | Biological replicates | OD\_timepoint

fig : plot curve of growth for each condition (ciblé et reel)

mes questions : est ce que replica bio, sont les meme entre stress et non stress ?voir avec Solène comment sont appliqué les stress (des debut ou apres un certain temps )=\> parler de ça en discussion =\> see annex for the media composition and more details

#table(
  columns: 3,
  align: (auto,auto,auto,),
  table.header([Col1], [Col2], [Col3],),
  table.hline(),
  [], [], [],
  [], [], [],
  [], [], [],
)
Tableau des conditions biologiques : - figure / tableau et explications des conditions biologiques

== microSPLiT protocol
<microsplit-protocol>
=== microSPLiT barcoding
<microsplit-barcoding>
#figure([
#box(image("../figures/protocol.png"))
], caption: figure.caption(
position: bottom, 
[
MicroSPLiT Protocol
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)
<fig-microsplit>


MicroSPLiT in-cell cDNA barcoding scheme. a, Bacterial cells are fixed overnight and permeabilized (Part 1, Steps 7–26) before the mRNA is preferentially polyadenylated (Part 1, Step 34). After mRNA enrichment, cells may contain both polyadenylated and non-polyadenylated mRNA. b, Cells are distributed into the first barcoding plate, and the mRNA is reverse transcribed by using a mixture of poly-dT and random hexamer primers carrying a barcode (barcode 1, BC1) and a 5’ phosphate for future ligation at their 5’ end (Part 1, Step 41). After the barcoding reaction, cells are pooled together and split again into the second barcoded plate (Part 1, Steps 43–52). c, Ligation adds a 5’ phosphorylated barcode 2 (BC2) to BC1 with a linker strand. A blocking solution is then added to each of the wells of the second plate, preventing any unreacted BC2 from future ligation (Part 1, Step 56). Cells are pooled and split into the third and final barcoded plate (Part 1, Steps 58–61). d, A second ligation step adds barcode 3 (BC3) with another linker strand. BC3 also contains a 5’ biotin, a primer binding site and a unique molecular identifier (UMI). A blocking solution for the R3 linker is added to each of the wells in the plate before the final pooling of cells (Part 1, Step 64). This results in uniquely barcoded cells that can be distributed in aliquots into sub-libraries and stored until future use or used immediately for library preparation. R1, round 1; R2, round 2; R3, round 3.

- Microbial split-pool ligation transcriptomics (microSPLiT) protocol was performed as described in @gaisser2024 and @kuchina2021. It’s a high-throughput single-cell RNA sequencing method for bacteria can profile transcriptional states in hundreds of thousands of bacteria in a single experiment without specialized equipment.

As bacterial samples are fixed and permeabilized before barcoding, they can be collected and stored ahead of time.

Contrary to other single-cell RNA sequencing methods, microSPLiT does not require the isolation of individual cells.

Instead of lysing bacteria and releasing the transcripts from each cell into a barcoding reaction vessel, in microSPLiT, each cell is the vessel enclosing its own transcripts. microSPLiT preserves the intact bacterial cell as a reaction compartment, enabling in situ barcoding of intracellular RNA.

==== Fixation and permeabilization
<fixation-and-permeabilization>
===== Fixation
<fixation>
For each sample, Bacterial cells are collected in bulk and fixed using formaldehyde. Fixation has two essential roles: - It preserves the transcriptomic state of each cell at the time of sampling. - It creates covalent cross-links between RNA and intracellular proteins, thereby preventing RNA leakage during later processing. - Importantly, fixation must retain the physical integrity of the cells to ensure downstream single-cell analysis is possible.

===== Permeabilization
<permeabilization>
Permeabilization is the process of making the cell membrane permeable to allow the entry of molecules. In microSPLiT, the cell wall/cell membrane is permeabilized using a combination of detergents and enzymes. The aim is to permeabilize the cell envelope without disrupting cell structure. This allows external enzymes (e.g., poly(A) polymerase, reverse transcriptase, ligase) and oligonucleotides to enter the cell. A key balance must be achieved: sufficient permeabilization for enzyme access, but minimal structural damage to maintain single-cell resolution.

==== In-cell polyadenylation
<in-cell-polyadenylation>
After permeabilization, the transcripts in the fixed and permeabilized cells undergo in situ polyadenylation with the addition of Escherichia coli poly(A) polymerase (PAP) and ATP. This step enriches for mRNA in the total barcoded RNA pool because, under these conditions, PAP preferentially polyadenylates mRNA as opposed to rRNA.

During the first barcoding round, the fixed and permeabilized bacteria are distributed into a 96-well plate, where their transcripts are reverse transcribed into cDNA and labeled with the first well-specific barcode inside the cells. The cells are mixed and redistributed two more times into new 96-well plates, where the second and third barcodes are appended to the cDNA via in-cell ligation reactions. Finally, the cells are mixed and divided into aliquot sub-libraries, which can be stored until future use or prepared for sequencing with the addition of a fourth barcode. It takes 4 days to generate sequencing-ready libraries, including 1 day for collection and overnight fixation of samples.

microSPLiT barcoding

Instead of lysing bacteria and releasing the transcripts from each cell into a barcoding reaction vessel, in microSPLiT, each cell is the vessel enclosing its own transcripts. The procedure starts with the collection of cells in bulk and the fixation of the bacterial suspension with formaldehyde. It then proceeds with permeabilization by using sequential mild detergent and lysozyme treatments (Fig. 1a). Fixation is critical because it both preserves the cellular transcriptomic state and covalently cross-links the transcripts with the proteins inside the cells to prevent leakage after permeabilization.

The permeabilization step ensures that the externally supplied enzymes and oligonucleotides can access the RNA transcripts in the fixed intracellular milieu. While sufficient permeabilization is crucial to the efficiency of barcoding, it is also critical to preserve the physical integrity of the fixed cells to maintain the single-cell resolution of the method.

We emphasize that for a successful microSPLiT experiment, the cells, after permeabilization, must still exist as intact, individual units to permit several split and pool steps and hold together the cross-linked RNA. After permeabilization, the transcripts in the fixed and permeabilized cells undergo in situ polyadenylation with the addition of Escherichia coli poly(A) polymerase (PAP) and ATP. This step enriches for mRNA in the total barcoded RNA pool because, under these conditions, PAP preferentially polyadenylates mRNA as opposed to rRNA.

==== round barcoding
<round-barcoding>
on split les 18 samples , en 5 pour avoir des replicats techni donc on obtient 90 samples , permettra d’evaluer la variance technique

depot meme quantité de cellules dans chaque puit

===== round 1 barcoding
<round-1-barcoding>
les 90 samples sont reparties dans 90 puits distincts, chacun contenant un unique primer barcodé

In the next step, the cell suspension and each sample is distributed into a 96-well plate with uniquely barcoded primers in each well (Fig. 1b, round 1 (R1) reverse transcription (RT) working plate).

The mRNA is then converted to cDNA through in-cell RT with a mixture of barcoded poly(T) and random hexamer primers. cette etape de barcoding permet demarquer les cellules par condition

===== round 2 barcoding
<round-2-barcoding>
Cells are then pooled, washed and randomly redistributed into a new 96-well plate (round 2 (R2) ligation working plate) containing a second set of well-specific barcodes, which are appended to the first barcode on the cDNA through an in-cell ligation reaction (Fig. 1c). Because of the random cell distribution, there is a high chance that each well of the secondround plate will contain cells with a mixture of different first-round barcodes, creating diverse barcode combinations.

===== round 3 barcoding
<round-3-barcoding>
Cells are then pooled again, and a split-ligation-pool cycle is repeated for the second time. Cells are randomly distributed into a third 96-well plate (round 3 (R3) ligation working plate), which is loaded with barcoded oligonucleotides containing the third cell barcode annealed with a linker, a 10-base unique molecular identifier (UMI), a common PCR handle and a 5’ biotin molecule . en réalité ici 95 a la place de 96 pref separer

90#emph[96];95 = 820800 combinations de barcodes possibles =\> autant de cellules individuelles possibles

==== 
<section>
The pooled cells are washed, counted and divided into sub-libraries of variable sizes, which can be stored at -80 °C for \>=6 months before proceeding with sequencing library preparation. Dividing sub-libraries into aliquots has two main advantages. First, it allows fine control over the number of cells in the final sequencing libraries. The size of a sub-library can be chosen so that the number of cells that receive the same barcode combination by chance does not exceed the desired collision rate (Table 1). It also permits multiplexing several libraries, potentially even from different experiments, in a single sequencing run.

on choisi de sequencer la plus petite librairie, celle avec 3000 cellules afin d’avoir de maximiser la profondeur de sequençage et limiter le nombre de collision rate (0.34 quand 96#emph[96];96) et avoir un nombre suffisant de cellules pour avoir un signal

90 samples see annex for the plate with barcoded primers

=== microSPLiT sequencing library preparation
<microsplit-sequencing-library-preparation>
== other remarks
<other-remarks>
Bacterial culture (Biological conditions) Before to explain the microSPLiT the experimental methods are presented.

see the protocol of - deux conditions env : stress (low\_glucose\_low fer) / pas stress

- suivi temporel à 3 timepoints (peut etre discuter de ce point après car les mesures de do ont été fait a des temps précis et variation de DO entre les reelles et attendus )

- 3 replicats biologiques par condition

- 5 replicats techniques par replicat biologique

\=\> un nombre cellules visés de 3000 cellules au totale (voir annexes pour le choix cela )

\=\> voir annexe pour le plan de plaques

\=\> donc 90 echantillons differents

\=\> rep bio : voir si variation entre les populations =\> replica techn : permettre d’estimer si variance dans le nombre - utilisation du nombre d’UMIs dans le round1 pour estimer cela =\> voir resultats : pooler ensemble

dans l’hypothèse equilibre parfaite a noté que 33 cellules par conditions ce qui est tres faible =\> peut etre soumis à des variations individuels au qui pourrait rendre difficule evalutaion variation (tirage aléatoire, possible cellules avec faible activi ou inverse (desequilibre))

\=\> discussion : renvoie vers l’outils Shiny pour voir les conditions biologiques permet une visualisation plus rapide des donnees

== MicroSPLiT
<microsplit>
figure du protocole microsplit

- explication de la méthode microSPLiT @brettner2024@kuchina2021

- pour faire simple : differentes etape : fixation ; …;

- contrairement aux autres methodes, pas d’isolation individuelle des cellules

- 3 rounds de split-pool : le premier round on affecte les differentes conditions biologiques

- un 4eme round pour rajouté un UMI pour les differentes libraires de sequençage

\-envoie d’un pool de librairies à la plateforme de sequençage - sequençage de type NovaSeq™ X Plus par la plateforme GenoBIRD , et demulteplexé permet d’amelioré la qualité de sequençage

- renvoie vers le protocole de kuchina 2021 et bretner 2024 @brettner2024@kuchina2021 pour l’explication de la méthode en detail

- =\> discussion des limites de cette methode dans la partie discussion

- ce qui est importantes de comprendre c’est que cela repose sur un methodes mathematiques combinatoires mais pas methodes d’isolation en tant que tel .

- un grand nombre de cellules ont été barcodés : plusieurs dizaines ou 100aines de milliers mais seulemtnn pres de 3000 cellules ont été choisi pour le sequencçage : \
  voir annex pour le tableau de choix du nombre de cellules pour la librairies (est un compromis pour avoir suffissament de cellules potentielles mais pas trop pour ne pas avoir des librairies trop grandes qui pourrait entrainer un profondeur de sequençage pas suffisante ) =\> discussion sur le nombre de cellules choisi pour la librairie

== Librairies structures
<librairies-structures>
\-stucture de la librairie

- figure : Final librairies structure

- R1 contient la sequence d’interet

- R2 contient les barcodes

- polyA ou random\_hexamer

\-\> key point : dans chauque puit polyA et random

\=\> voir annex pour la structure complete avec TSO… =\> discussion sur TSO

\-STARsolo permet l’alignements des reads et lectures des barcodes - details de methode d’alignements partielles ou non … - tailles des reads que j’ai alignés

= Pipeline of the analysis
<pipeline-of-the-analysis>
- figure of the pipeline

Preprocesing sur le cluster genouest, la suite en local et sur le cluster aussi (peut etre) tous les scripts sont dispo sur differents depots githubs

- differentes etapes :

  - demultiplexage des index de librairies (réalisé par la plateforme de sequençage)

  - QC control des données avec Andrews, S. (2010). FastQC: ~A~Quality Control Tool for High Throughput Sequence Data \[Online\]. Available online at:~#link("http://www.bioinformatics.babraham.ac.uk/projects/fastqc/") et avec @ewels2016

  - trimming des données de seqeunçage avec Fastp @chen2018 et Cutadapt @martin2011

  - aligment des data sur le genome de reference de #emph[Pseudomonas brassicacearum] grace à STARsolo un amelioration de l’outils STAR pour les données single-cell \[#cite(<dobin2013>, form: "prose");\]@kaminow

  - different other tools existe comme pour alignement et lecture des barcodes SPLiTseq/ microSPliT comme Kallisto @bray2016, @sullivan2025 mais d’apres le benchmarking de le plus rapide, reproductible est starsolo @kuijpers2024

    Although single‐cell sequencing approaches have been developed for several molecular modalities, single‐cell transcriptome sequencing is the most prevalent and widely applied technique. SPLiT‐seq (split‐pool ligation‐based transcriptome sequencing) is one of these single‐cell transcriptome techniques that applies a unique combinatorial‐barcoding approach by splitting and pooling cells into multi‐well plates containing barcodes. This unique approach required the development of dedicated computational tools to preprocess the data and extract the count matrices. Here we compare eight bioinformatic pipelines (alevin‐fry splitp, LR‐splitpipe, SCSit, splitpipe, splitpipeline, SPLiTseq‐demultiplex, STARsolo and zUMI) that have been developed to process SPLiT‐seq data.

  - Metadata assignation to seura

\=\> d’autres étapes ou autres outils pourrait etre ajouter dans le pipeline (voir la partie disccussion)

- differnetes steps of the pipeline
-
- assignation des metadonnées (utilisation d’un genome de reference =\> discussion )

=== Demultiplexing
<demultiplexing>
=== 
<section-1>
=== 
<section-2>
=== STARsolo
<starsolo>
500 go

attention je vais devoir dire les version des outils utilisés

github,

differentes etapes mais pas detailler tous ici renvoie vers le codes commenter et les readme pour comprendre en details telechangement des fichiers, decompression .zip ; 4 fichier , analyse indenpendant de la qualité des 4 librair, trimming

parallelisation pour ganger du temps …

===== Final effective command line of STARsolo:
<final-effective-command-line-of-starsolo>
#strong[STAR] \
#strong[–runThreadN] 64 \
#strong[–genomeDir] /path/to/genome\_index \
#strong[–readFilesIn] \
\/path/to/input/merged\_trimmed-R1.fastq.gz \
\/path/to/input/merged\_trimmed-R2.fastq.gz \
#strong[–readFilesCommand] gunzip -c \
#strong[–outFileNamePrefix] /path/to/output/starsolo\_output/ \
#strong[–outSAMtype] BAM Unsorted \
#strong[–outFilterScoreMinOverLread] 0 \
#strong[–outFilterMatchNmin] 50 \
#strong[–outFilterMatchNminOverLread] 0 \
#strong[–alignSJoverhangMin] 1000 \
#strong[–alignSJDBoverhangMin] 1000 \
#strong[–soloType] CB\_UMI\_Complex \
#strong[–soloCBwhitelist] \
\/path/to/barcodes/barcode\_round3.txt \
\/path/to/barcodes/barcode\_round2.txt \
\/path/to/barcodes/barcode\_round1.txt \
#strong[–soloFeatures] Gene GeneFull \
#strong[–soloUMIdedup] 1MM\_All \
#strong[–soloCBmatchWLtype] 1MM \
#strong[–soloCBposition] 0\_10\_0\_17 0\_48\_0\_55 0\_78\_0\_85 \
#strong[–soloUMIposition] 0\_0\_0\_9 \
#strong[–soloMultiMappers] Uniform

=== STARsolo Parameters Explanation
<starsolo-parameters-explanation>
This section details the key parameters used in our STARsolo analysis and their significance:

==== General STAR Parameters
<general-star-parameters>
- `--runThreadN 64` : Use of 64 threads for parallel alignment
- `--genomeDir` : Path to the reference genome index
- `--readFilesIn` : Input FASTQ files (R1 and R2)
- `--readFilesCommand gunzip -c` : Command to decompress FASTQ.gz files
- `--outFileNamePrefix` : Prefix for output files
- `--outSAMtype BAM Unsorted` : Unsorted BAM output format

==== Filtering Parameters
<filtering-parameters>
- `--outFilterScoreMinOverLread 0` : Minimum filtering score relative to read length
- `--outFilterMatchNmin 50` : Minimum number of matching bases for a valid alignment
- `--outFilterMatchNminOverLread 0` : Minimum match ratio relative to read length
- `--alignSJoverhangMin 1000` and `--alignSJDBoverhangMin 1000` : Strict parameters for splice junction detection

==== STARsolo-specific Parameters
<starsolo-specific-parameters>
- `--soloType CB_UMI_Complex` : Analysis type for cell barcodes (CB) and complex UMIs
- `--soloCBwhitelist` : List of valid cell barcodes for the three barcoding rounds
- `--soloFeatures Gene GeneFull` : Analysis of features at both gene and full transcript levels
- `--soloUMIdedup 1MM_All` : UMI deduplication with one mutation tolerance
- `--soloCBmatchWLtype 1MM` : Cell barcode matching with one mutation tolerance
- `--soloCBposition` : Cell barcode positions in reads (3 rounds)
  - Round 1: 0\_10\_0\_17
  - Round 2: 0\_48\_0\_55
  - Round 3: 0\_78\_0\_85
- `--soloUMIposition 0_0_0_9` : UMI position in reads
- `--soloMultiMappers Uniform` : Uniform distribution of multi-mapped reads

These parameters were chosen to optimize single-cell detection while maintaining high alignment quality and accounting for the complexity of our three-round barcoding protocol.
