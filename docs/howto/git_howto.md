# Git

General Git notes, mostly for our own benefit

- - -

add remote tracking branch (`upstream/foo`) to current branch

    git branch -u upstream/foo

add remote tracking branch (`upstream/foo`) to another branch `foo`

    git branch -u upstream/foo foo

- - -

list remotes

    git remote -v

- - -

list all details for a specific remote

    git remote show <remote_name>

- - -

list all branches (including remotes)

    git branch -a

- - -

git detailed history for a single file

    git log --follow <filename>

- - -

switch to a branch

    git checkout <branch_name>

- - -

create a new branch from the current branch and switch to it

    git checkout -b <branch_name>

- - -

create a new branch from another branch and switch to it

    git checkout -b <branch_name> <another_branch>

- - -

push a local branch to a remote (`origin`)

    git push -u origin <branch_name>

- - -

delete a local branch (after it has been merged)

    git branch -d <branch_name>

- - -

delete a local branch (before it has been merged)

    git branch -D <branch_name>

- - -

delete a remote branch

    git push origin --delete <branch_name>

- - -

undo, edit and redo the last commit

    git reset --soft 'HEAD^'
    [edit the commit]
    git add [files to be staged]
    git commit -c ORIG_HEAD

or if you do not need to edit the original commit message, use

    git commit -C ORIG_HEAD

- - -

amend a commit

    git commit --amend

- - -

undo the last commit (so it will be gone)

    git reset --hard 'HEAD^'

or

    git reset --hard HEAD~1

- - -

abort a merge

    git merge --abort

- - -

diff two branches

    git diff <branch_1>..<branch_2>

will produce the diff between the tips of the two branches. To find the diff from their common ancestor, use three dots instead of two:

    git diff <branch_1>...<branch_2>

# ctags and cscope

create ctags database

    ctags -R -f .tags

- - -

cscope gather files and create database

    find . -name '*[ch]' -not -path "*git*" -type f > cscope.files
    cscope -b -q