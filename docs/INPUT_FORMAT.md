# Trait types and units

Keep the updated observations on the first worksheet. Put metadata labels in the Variety
column (not Entry or the Candidate/Check Type column), below the observations:

| Variety | Plant vigor | Fruit weight | Fruits/plant | Fruit set |
| --- | --- | --- | --- | --- |
| DIRECTION | high | high | high | high |
| WEIGHT | 3 | 5 | 4 | 3 |
| TYPE | score | quantitative | quantitative | quantitative |
| UNIT | score | kg | number | percent |

- TYPE is optional: score, quantitative, continuous or continues (case insensitive).
- UNIT is optional: score, gr, kg, number, cm, percent, brix, or your measurement unit.
- UNIT alone can declare the type: score means scoring; kg/cm/etc. mean quantitative.
- A score can contain fractions such as 3.5, but should remain within 1-5.
- Percent means 0-100. Out-of-range values are flagged for breeder review.
- Number may represent averaged counts, so decimals are not automatically rejected.
- Leave units blank when unknown. Quantitative/continuous in UNIT also works for
  simple workbooks, but TYPE is the preferred place for these labels.
- Do not mix a score TYPE with a kg/cm/etc. UNIT. Use one TYPE and one UNIT row.
- With neither row declared, the previous inference remains and is labelled inferred.
- Entry/ID is a plot identifier, not a trait. Variety can repeat across replications.
  Repeated Entry/ID values generate a warning; they are not automatically removed.
- 2500 means 2500; 2.5 and 2,5 both mean 2.5. A dot/comma never means thousands.
  Mixed/repeated separators become missing numeric values and should be corrected.
  Excel may convert text before import: check the actual stored numeric value.

After upload, the summary reports scoring and quantitative trait counts.
The upload summary names the file/worksheet and reports detected observation-ID
columns and TYPE/UNIT row counts. The preview shows that sheet unchanged; it does
not generate an ID column or invent missing metadata. Save Excel edits and reupload
the updated file. The first worksheet is read automatically.
Metadata is included in full analysis workbook exports.
Declaring a type does not by itself change the experimental-design model or add
an ordinal-model implementation. Units are preserved, not automatically converted.
