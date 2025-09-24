---
name: 📖 User Story
about: Capture analytical requirements as user stories with clear acceptance criteria using Gherkin syntax.
title: ""
labels: "user-story, needs-triage"
assignees: "ASLM-Fabebe"
---

**As a** [researcher/clinician/lab technician/data analyst]  
**I need** [specific functionality]  
**So that** [business/research value]  

### Details and Assumptions  
* [Document known constraints, dependencies, data requirements]
* [Specify any regulatory or compliance considerations]
* [Note integration points with existing AMR/AMC/AMU workflows]

### Acceptance Criteria  
Primary scenario
```gherkin
Given [initial system state and user context]
When [user takes specific action]
Then [observable outcome occurs]
And [additional verification criteria]
```

Edge cases or alternative flows
Scenario: [Alternative scenario name]
```gherkin
Given [different starting conditions]
When [user action]
Then [expected behaviour]
```

### Data Requirements
_If applicable:_
- [ ] Requires sample data for testing
- [ ] Works with existing test datasets
- [ ] Needs specific data format/structure
- [ ] Has regulatory data handling requirements

### Regulatory Context
- [ ] Supports GLASS reporting
- [ ] Aligns with WHO surveillance guidelines
- [ ] Required for CLSI/EUCAST compliance  
- [ ] Other: ________________
