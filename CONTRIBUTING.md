# Contributing to the G4CMP Package

Last updated 25 September 2024, Michael Kelsey


The "Geant4 Condensed Matter Physics" (G4CMP) package is open source (see
[LICENSE](LICENSE)), but we are limiting code contributions to registered
developers from particle physics or related experimental collaborations.  If
you'd like to contribute, please send an e-mail to the package owner,
currently Michael Kelsey (Texas A&M) <kelsey AT slac.stanford.edu>, to be
added to the Contributors list for G4CMP.


## Developer Consortium

The [G4CMP Developers
Consortium](https://confluence.slac.stanford.edu/display/G4CMP/G4CMP+Developer%27s+Consortium)
provides coordination and communication among G4CMP developers and users.
Joining the Consortium also provides developer access to the G4CMP
repository to create feature branches, issues, and pull requests, as
described below.


## Reporting Problems

Until recently, G4CMP has been used primarily (and almost exclusively) by the
SuperCDMS Collaboration, members of which originally developed G4CMP.  As a
result, we have been using an internal issue tracking system in JIRA for
problems, requests, and development work.  We have made this tracking system
public on the [SLAC JIRA server](https://jira.slac.stanford.edu/browse/G4CMP);
G4CMP users may apply for a Crowd Account for access.


## Feature Tracking

The [G4CMP JIRA project](https://jira.slac.stanford.edu/browse/G4CMP) is
used to track feature development as well as problem reporting.  Please
create a JIRA ticket _before_ you begin developing a new feature; assign
yourself to the ticket.  When you start work, please modify the ticket state
from "Open" to "In Progress," and provide descriptive comments as you go
along.  The JIRA ticket number will be used (see below) in naming the
feature branch containinhg your work, and will be used to tag your
development when they are merged onto **develop** branch.


## Contribution Workflow

We are using a simplified version of the ["Git Flow"
workflow](https://nvie.com/posts/a-successful-git-branching-model/).  There
is a single line of development, on the **develop** branch, and a single
line of tagged "releases" (see below) on the **master** branch.  We are not
currently planning to rename the latter, to avoid disrupting existing clones
of the repository out in the world.


## External Forks

We would prefer that external forks not be used for submitting software
changes.  If you have created a fork, please create a JIRA ticket and a
feature branch within the repository, as discussed below.  We can then use a
PR to merge your fork onto the new feature branch, before completing the
process with our normal workflow.


## New Features

Contributors should make a named feature branch in this repository (rather
than a fork), branching from the HEAD of **develop** branch:
```
  git checkout develop
  git pull
  git checkout -b <branch-name>
```

Contributors should push their code changes directly back to the repository,
on the feature branch created above.  Do **not** commit or push changes
directly to the develop branch, and do **not** use the GitHub Web interface
to directly modify files.

```
  [ edit, edit, edit ]
  git commit -m "Explanatory message" ...
  git push -u
```

Our style guide for code contributions is 
located on this Confluence page: ["G4CMP Style Guide"](
https://confluence.slac.stanford.edu/spaces/G4CMP/pages/567218502/G4CMP+Style+Guide)
.


## ChangeHistory

The [ChangeHistory](ChangeHistory) file is used to track new features, bug
fixes, and other changes to repository.  It contains single-line dated
entries for each feature branch, named using either the [JIRA
ticket](https://jira.slac.stanford.edu/browse/G4CMP/) or the GitHub issue/PR
number associated with the feature.  

The entries are presented in reverse order (the most recent changes are at
the top of the file), and each **master**-branch release is identified with
a date and tag string A section of the file, before and after a release tag,
might look like this:

```
2024-05-02  G4CMP-378 : Correct expression for phonon-QP scattering energy.
2024-05-02  G4CMP-344 : Improve KaplanQP memory usage with mutable buffers.

2024-04-29  g4cmp-V08-09-00 Support manual setting of voltage bias for sampling.
2024-04-20  G4CMP-403 : Allow clients to pass voltage bias through HitMerging.

2024-04-19  g4cmp-V08-08-00 Fix major global/local bug, compiler warnings.
2024-04-18  G4CMP-402 : Reduce overhead generaing track positions along steps.
```

When a feature branch is ready for merging, the contributor should add a
dated line to the top of the [ChangeHistory](ChangeHistory) file, with the
branch name and a brief explanation of its purpose.  The contributor should
then submit a Pull Request on GitHub, or send e-mail directly to the package
owner (above), requesting that it be merged onto **develop** branch.


## Large Scale Projects

If you have started, or are participating in, a large scale project, we
strongly recommend using JIRA's "sub-task" system to track the many
individual development tasks within the project.  Start by creating an
overall JIRA ticket, and create a feature branch with that top-level name,
and update the ChangeHistory file to support tracking the whole project.

At the top of the ChangeHistory file, create a block of text as follows:

```
[ Modifications included on branch XYZ, put subtask tags below it: ]
2024-13-35  XYZ       : Create a big new feature in G4CMP.         

[ Main line "develop" ChangeHistory; when merging, put develop updates here: ]
```

followed by the most recent "develop" record.  Commit and push that change
on the top-level feature branch back to the repository.

For individual items within the project, create feature branches from the
top-level feature branch, rather than directly from **develop**.  Follow the
same workflow model as below, but when the individual task is completed,
merge it back onto the top-level feature branch with a tag.  In the
ChangeHistory file, insert the records for each task below the top-level
line above (but keeping the tasks in reverse date order).  Update the
datestamp on the top-level line to match the merge, for example:

```
[ Modifications included on branch XYZ, put subtask tags below it: ]
2025-02-01  XYZ       : Create a big new feature in G4CMP.
2025-01-01  G4CMP-997 : Generate a New Year's energy deposit in diamond.
2024-13-35  G4CMP-995 : Add function to compute phonon mass in diamond.

[ Main line "develop" ChangeHistory; when merging, put develop updates here: ]
2024-08-22  G4CMP-423 : Avoid wrong volume assignment in G4CMPHitMerging.
```

When the large scale project is completed, then it would be merged onto
**develop** branch as described below.  The ChangeHistory file can be edited
very easy to simply remove the bracketed instructions above, leaving what
looks like a normal development sequence:

```
2025-02-01  XYZ       : Create a big new feature in G4CMP.
2025-01-01  G4CMP-997 : Generate a New Year's energy deposit in diamond.
2024-13-35  G4CMP-995 : Add function to compute phonon mass in diamond.
2024-08-22  G4CMP-423 : Avoid wrong volume assignment in G4CMPHitMerging.
```


## Merging Onto Develop

When a feature branch is merged onto **develop** branch, the package owner
will tag the merge commit with the name of the feature branch, and (first)
delete the feature branch from the repository.
```
  git checkout <branch-name>
  git pull
  git checkout develop
  git pull
  git merge --no-ff -m "Explanatory message" <branch-name>
  [ resolve conflicts, git merge --continue, etc. ]
  git branch -d <branch-name>
  git push origin :<branch-name>
  git tag -a -m "Explanatory message" <branch-name>
  git push ; git push --tag
```


## Pseudo-Release Tags

We have not been using GitHub's "Release" capability.  Instead, we have only
been merging the **develop** branch onto **master**, and applying a
[semantic verioning](https://semver.org) tag there, of the form
**g4cmp-Vxx-yy-zz**.  This is done exactly as with Merging onto Develop
(above), but without deleting the **develop** branch.

Before merging, the [ChangeHistory](ChangeHistory) file should be updated
with the new release tag and a brief explanation.  The individual feature
branches include will then be listed chronoligally below that line.
