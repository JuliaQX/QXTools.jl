# Contributing to QuantEx

We welcome contributions and here we lay out some guidelines which should be followed to make the process more streamlined for all involved.

## Contribution process

Commits should not be pushed directly to the master branch but should instead be merged from feature branches following via pull/merge requests. To track 
tasks, features, bugs and enhancements should have a corresponding issue which explains the motivation and logic used in the committed code/documentation.
The steps in the full process from creating an issue to merging are:

1. Create issue with quick description of feature, bug, enhancement etc.. This will be assigned an issue number
2. Create a branch with a name that starts with the issue number and gives a concise description of the issue
3. Make the necessary changes to the branch. Commit messages should follow imperative style ("Fix bug" vs "Fixed bug"). Further guidelines for commit messages [here](https://gist.github.com/robertpainsi/b632364184e70900af4ab688decf6f53)
4. Create a merge/pull request requesting to merge the changes into the master branch and select the appropraite merge/pull request template, prefix the issue name with `WIP:` to indicate that it is a work in progress
5. The merge/pull request template has a number of check list items which should be satisfied before the request can be merged. These include:
- All discussions are resolved
- New code is covered by appropriate tests
- Tests are passing locally and on CI
- The documentation is consistent with changes
- Any code that was copied from other sources has the paper/url in a comment and is compatible with the MIT licence
- Notebooks/examples not covered by unittests have been tested and updated as required
- The feature branch is up-to-date with the master branch (rebase if behind)
- Incremented the version string in Project.toml file
6. Once all the items on the checklist have been addressed the `WIP:` prefix should be removed and the merge/pull request should be reviewed by one of the developement team and merged if accepted or comments added explaining rationale if not