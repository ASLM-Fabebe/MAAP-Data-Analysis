---
name: 🐛 Bug Report
about: Report an error or unexpected behaviour in AMR/AMC/AMU analysis. Replace placeholders with actual values and use the Given/When/Then format.
title: "[Bug]: "
labels: "bug, needs-triage"
assignees: "ASLM-Fabebe"
---

**Analysis Module:**  
_Select one:_ AMR / AMC / AMU / Setup / Data Processing

**Bug Description**  
_What went wrong? Be specific about the expected vs actual behaviour._

**Steps to Reproduce**  
_Provide minimal reproducible steps:_
1. Navigate to test-data folder and add file: `[filename]`
2. Run command: `library(shiny); runApp('[script_name]')`
3. Select options: `[specific selections]`
4. Observe error at step: `[which step]`

**Environment Details**  
- **OS:** [Windows 10/11, macOS, Ubuntu, etc.]
- **R Version:** [e.g., R 4.4.3]
- **RStudio Version:** [if applicable]
- **AMR Package Version:** [run `packageVersion("AMR")`]
- **Data File:** [e.g., `AMR_kenya_data.xlsx`, `AMC_test_data.xlsx`]
- **File Size/Records:** [approximate number of rows]

**Sample Data Available**  
- [ ] I can provide a sanitized sample dataset that reproduces this issue
- [ ] This issue occurs with the provided test data
- [ ] This issue requires specific data that cannot be shared

**Regulatory/Compliance Impact**  
- [ ] This affects GLASS reporting requirements
- [ ] This impacts surveillance data quality
- [ ] This affects regulatory compliance (specify): ________________
- [ ] No regulatory impact

**Acceptance Criteria** 
Example - Replace with your specific scenario
```gherkin
Given I have AMR data loaded in the correct format
When I run the analysis with standard parameters
Then the analysis completes without errors
And I receive the expected output files in the results folder
And all visualizations render correctly
```


**Error Details**  
_Console output, error messages, or screenshots:_
[Paste error messages here]


**Additional Context**  
_Any other relevant information, screenshots, or workarounds attempted._
