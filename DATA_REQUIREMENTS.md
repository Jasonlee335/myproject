# Data requirements

No participant-level inputs are included. Obtain data from the official providers:

- NHANES: https://wwwn.cdc.gov/nchs/nhanes/
- KNHANES: https://knhanes.kdca.go.kr/knhanes/
- CHARLS: https://charls.charlsdata.com/
- Harmonized CHARLS: https://g2aging.org/ and the source-data access conditions.

The archived analysis uses NHANES 2005-2010, 2013-2014, and 2017-March 2020;
KNHANES 2008-2011 and 2024; and CHARLS follow-up from 2011 through 2020.
Use the exact file/variable definitions in `01_data_cleaning.R`; do not substitute
survey years or alter variable coding to mimic unavailable fields.

## Expected local source layout

`DXA_DATA_ROOT` must contain the existing filenames referenced by Step 01:

```text
nhanes_osteoporosis_screening/raw/
KNHANES/raw/
CHARLS.csv
charls_official/raw_extracted/2011/
charls_official/raw_extracted/2013/
charls_official/raw_extracted/2015/
charls原版数据/Harmonized_CHARLS_D/H_CHARLS_D_Data.dta
charls原版数据/2018/CHARLS2018r/
charls原版数据/2020/CHARLS2020r/
```

The non-ASCII directory name is retained to match the archived input path. If you
rename folders, deliberately update the corresponding local path definitions.
Baseline modules support the FRAX metadata audit; weights and exit modules support
the longitudinal and mortality checks.

## Consolidated CHARLS input

`CHARLS.csv` is a pre-existing author-prepared longitudinal table, not an official
release filename. Create or obtain it lawfully. Upstream consolidation code is not
in the current archive; this is not a complete raw-to-longitudinal harmonization
package. Columns consumed by Step 01 are:

| Column | Meaning or consolidated coding |
| --- | --- |
| ID | Participant identifier as text |
| householdID | Household identifier as text |
| communityID | Community identifier as text |
| wave | Survey wave 1 through 5 |
| gender | 0 female, 1 male |
| age | Years |
| mweight | Measured weight, kg |
| mheight | Measured height, cm |
| bmi | Consolidated BMI; calculated BMI is a fallback |
| hip | Hip-fracture indicator, 0/1 |
| fall_down | Falls indicator, 0/1 |
| province | Province identifier or label |
| rural | Residence indicator, 0/1 |
| stroke | Stroke indicator, 0/1 |
| arthre | Arthritis/rheumatism indicator, 0/1 |
| smoken | Current-smoking indicator, 0/1 |
| drinkev | Ever-drinking indicator, 0/1 |
| exercise | Exercise indicator, 0/1 |

These are consolidated input codes, not a claim that all raw questionnaires use
the same coding. Preserve missing values and questionnaire meanings. One row per
participant-wave is required; duplicate keys stop execution.

The optional `charls_official/audit_original/key_variable_validation.csv` supports
a prior source-count check when available. Archived input hashes identify the
clean inputs but do not provide access to them.
