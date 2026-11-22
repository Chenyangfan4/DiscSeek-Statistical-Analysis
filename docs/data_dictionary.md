# Input data dictionary

No patient-level data or example patient records are distributed. Input CSV headers must match the exact names below. Standard `read.csv()` converts spaces and punctuation to periods in R column names; downstream scripts use those converted names.

## External prospective validation: `data/external/total_ivd.csv`

Each row represents one cervical intervertebral disc level of a selected MRI examination. One selected examination contributes multiple disc-level rows (C2/C3–C7/T1). Age and sex are demographic attributes, not image-derived measurements. IVDType is a disc-level label, not a measurement.

| Source CSV column | Expected analysis type | Values / levels | Role or definition |
|---|---|---|---|
| ScanSite | Factor | QQHE; DT; ZJ | Acquisition hospital (see study hospital labels) |
| ScanMachine | Factor | Philips; GE; Siemens; UIH | MRI scanner vendor |
| ScanType | Factor | Manual; AI; Interactive | MRI scan planning approach |
| ImageID | Identifier |  | Unique identifier of the selected MRI examination |
| PatientName | Identifier |  | Confidential source field; excluded from processed outputs |
| PatientID | Identifier |  | Source patient identifier; excluded from processed outputs |
| Age | Numeric |  | Patient age at examination, years |
| Sex | Factor | F; M | Patient sex as recorded in clinical data |
| Max.scan.angle.offset | Numeric |  | Maximum scan angle offset within an examination (degrees) |
| IVDType | Factor | C2C3; C3C4; C4C5; C5C6; C6C7; C7T1 | Cervical disc level represented by the current row |
| Grade.of.IVD.protrusion | Ordered categorical |  | See Supplementary Table S7 |
| Grade.of.left.foraminal.stenosis | Ordered categorical |  | See Supplementary Table S7 |
| Grade.of.right.foraminal.stenosis | Ordered categorical |  | See Supplementary Table S7 |
| Grade.of.spinal.stenosis | Ordered categorical |  | See Supplementary Table S7 |
| Height.of.IVD..cm. | Numeric |  | See Supplementary Table S8; cm |
| Height.of.left.lateral.recess..mm. | Numeric |  | See Supplementary Table S8; mm |
| Height.of.right.lateral.recess..mm. | Numeric |  | See Supplementary Table S8; mm |
| Morphology.symmetry.of.bilateral.lateral.recesses | Categorical |  | See Supplementary Table S7 |
| Sagittal.diameter.of.IVD..cm. | Numeric |  | See Supplementary Table S8; cm |
| Sagittal.diameter.of.VB..cm. | Numeric |  | See Supplementary Table S8; cm |
| Sagittal.diameter.of.spinal.canal..cm. | Numeric |  | See Supplementary Table S8; cm |
| Scan.angle.offset | Numeric |  | Disc-level axial scan angle offset (degrees) |
| Torg.ratio | Numeric |  | See Supplementary Table S8; dimensionless |
| Transverse.diameter.of.IVD..cm. | Numeric |  | See Supplementary Table S8; cm |
| Transverse.diameter.of.spinal.canal..cm. | Numeric |  | See Supplementary Table S8; cm |
| Radiology.diagnosis | Factor | ACS; CSM; CSR | Radiological disease-subtype classification |
| Clinical.diagnosis | Factor | ACS; CSM; CSR | Clinical disease-subtype classification |
| Imaging.quality.score | Ordinal | Observed 2; 3; 4; 5 | Image quality rated on a five-point scale; observed levels are 2–5 |

The imaging-quality instrument allows scores 1–5; the supplied observed levels are 2–5. Anatomical grading and measurement rules are specified in Supplementary Tables S7 and S8 and are not redefined here.

### Data protection and grouping

The processed analysis file omits direct personal identifiers. `ImageID` identifies the selected examination and is used for held-out partitioning. Four overlapping groups are constructed only for descriptive summaries; the inferential analyses preserve nonduplicated examinations.

## Internal technical validation: `data/internal/internal_planning_error.csv`

Each record is one MRI examination from the independent internal technical testing cohort. The following fields were described in the original `internal_analysis.R` script.

| Column | Definition | Unit |
|---|---|---|
| ID | Internal test examination identifier | — |
| SAG_dirET | Sagittal inter-technologist scan angle difference | degrees |
| SAG_cenEA | Sagittal DiscSeek-to-expert FoV center offset | mm |
| SAG_dirEA | Sagittal DiscSeek-to-expert scan angle offset | degrees |
| COR_dirET | Coronal inter-technologist scan angle difference | degrees |
| COR_cenEA | Coronal DiscSeek-to-expert FoV center offset | mm |
| COR_dirEA | Coronal DiscSeek-to-expert scan angle offset | degrees |
| AXI_dirET | Axial inter-technologist scan angle difference | degrees |
| AXI_cenEA | Axial DiscSeek-to-expert FoV center offset | mm |
| AXI_dirEA | Axial DiscSeek-to-expert scan angle offset | degrees |
| SAG_DISC | Visible cervical disc count in best-viewed sagittal image | count |
| COR_DISC | Visible cervical disc count in best-viewed coronal image | count |
| AXI_DISC | Visible cervical disc count in scout MRI | count |
| Qualified.Mask | Technologist assessment of whether the disc segmentation mask was sufficient for scan prescription | category |
| Quality.Factor | Calculated score for segmentation mask qualification | score |

## V-Net training logs: `data/training/loss_metric_curves.csv`

Each row records a training or validation metric at a particular optimization iteration or epoch.

| Column | Description |
|---|---|
| Epoch | Training epoch index |
| Iteration | Training or validation iteration index |
| Loss | Loss value for the specified network and dataset split |
| Dice_ID_Rate | Segmentation Dice or identification rate, depending on the network |
| Dataset | Training/validation split label used by the source logging code |
| Network | Segmentation or identification V-Net label used by the source logging code |

