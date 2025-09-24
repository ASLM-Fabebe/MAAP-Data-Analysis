---
name: 📚 Documentation Improvement
about: Suggest updates or corrections to documentation. Define clear criteria for when documentation is considered complete.
title: "[Docs]: "
labels: "documentation, needs-triage"
assignees: "ASLM-Fabebe"
---

**Documentation Type:**  
_Select one:_
- [ ] **README** - Main repository instructions
- [ ] **Installation Guide** - R/RStudio/Package setup
- [ ] **AMR Analysis Guide** - Set A documentation  
- [ ] **AMC Analysis Guide** - Set B documentation
- [ ] **AMU Analysis Guide** - Set C documentation
- [ ] **Data Format Guide** - Input file specifications
- [ ] **Troubleshooting Guide** - Common issues and solutions
- [ ] **Code Comments** - In-script documentation
- [ ] **Output Interpretation** - Results explanation
- [ ] **API/Function Documentation** - Technical reference

**Current Location**  
_Specify file path, URL, or section:_

**Issue Description**  
_What is unclear, missing, outdated, or incorrect?_

**Suggested Improvement**  
_How should it be updated? Include specific text changes if possible._

**Missing Gherkin Examples**  
- [ ] Documentation lacks clear acceptance criteria examples
- [ ] Need Gherkin scenarios for common use cases
- [ ] Should include Given/When/Then examples for troubleshooting

**Target Audience Impact**  
_Who struggles with the current documentation?_
- [ ] **New Users** - First-time setup and basic usage
- [ ] **R Developers** - Technical implementation details
- [ ] **Clinical Researchers** - Interpretation and application
- [ ] **Laboratory Personnel** - Data preparation and workflows
- [ ] **Data Analysts** - Advanced analysis techniques
- [ ] **System Administrators** - Installation and maintenance

**Acceptance Criteria**  
Documentation completeness
```gherkin
Given the updated documentation
When a [target user type] follows the instructions
Then they successfully complete [specific task]
And they understand [key concepts]
And they can troubleshoot common issues independently
```
Clarity verification
```gherkin
Given a new user reads the documentation
When they encounter [specific scenario]
Then the information is clear and actionable
And they don't need to consult external resources
```

**Validation Method**  
_How can we verify the documentation improvement works?_
- [ ] Test with a new user following the updated instructions
- [ ] Validate against common support questions
- [ ] Check completeness against similar projects
- [ ] Review with domain experts (clinical/laboratory staff)

**Related Issues/Discussions**  
_Link any related issues, discussions, or support requests._
