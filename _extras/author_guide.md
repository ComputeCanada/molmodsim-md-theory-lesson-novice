---
layout: page
title: "Authors Guide"
permalink: /author-guide/
questions:
- "How can I create lessons using this style?"
- "How do I update the style of my lesson?"
- "How can we add mathematical equations into our lessons?"
- "How can we add Flow diagrams as code within the lessons?"
objectives:
- "Give a quick overview on how this style is used."
- "Demonstrate the use of MathJax equations inside a lesson."
- "Demonstrate the use of Mermaid charts inside a lesson."
keypoints:
- "Use `$$ LaTeX math $$` within a text."
- "Use `$ LaTeX math $` followed by `{: .math-center}` or `{: .math-left}` for math-blocks."
- 'Use a `<div class="mermaid"></div>` to add the code for a graph.'
---

* Table of Contents
{:toc}

# Preface

[This style]({{ site.cc_style_repo }}) is based on based on The Carpentries style and has some 
modifications for use by [Compute Canada]({{ site.cc_site }}) and [ACENET]({{ site.an_site }}). 
Creating and editing lessons from such styles is explained in the 
[Carpentries Lesson Example]({{ site.example_site}}).

# Creating a Lesson With This Style

Here is a quick rundown on how to create a lesson based on this style.

## Initializing a lesson

~~~
$ git clone {{ site.cc_style_repo }}.git  my_lesson
$ cd my_lesson
$ python3 bin/lesson_initialize.py
$ git remote rename origin template
$ git remote add origin  <GITHUB_URL_FOR_YOUR_LESSON>
$ git commit --all  -m "initialize lesson"
$ git push -u  origin  gh-pages
~~~
{: .language-bash }

## Configuring the Lesson

Edit the file `_config.yml` to set these values:
* `carpentry`:  This changes the branding of the lesson, e.g. `carpentry: "cc"`
   for Compute Canada.
* `title`: Sets the lesson title.
* `life_cycle`: This indicates how far along the development of the lesson has come.
* `email`: The email address of the primary contact (probably an email list).

Also don't forget to edit these files:

* `AUTHORS`
* `CITATION`
* `index.md`
* `README.md`
* `reference.md`
* `setup.md`
* add your episodes (chapters) under `_episodes/*`


## Updating the style

> ## Make sure to resolve conflicts!
> 
> Merging changes from the template can result in merge conflicts.
> 
> Proficiency in resolving conflicts is recommended, as is having Jekyll set up to test the page
> locally.  But as we are working with git repositories, we can always go back to an earlier version
> of the history.
{: .callout }

~~~
# check if you have a remote "template" and whether it points to {{ site.cc_style_repo }}.git :
$ git remote -v
origin	git@github.com:ComputeCanada/my_lesson.git (fetch)
origin	git@github.com:ComputeCanada/my_lesson.git (push)
template	{{ site.cc_style_repo }}.git (fetch)
template	{{ site.cc_style_repo }}.git (push)

# if you don't already have the remote "template", add it:
$ git remote  add  template  {{ site.cc_style_repo }}.git

# now "fetch" the latest commits from "template":
$ git fetch template

# merge the changes:
$ git merge template/gh-pages

# Now resolve the conflicts and test if the page still works.
~~~
{: .language-bash }

# Using Syntax Extensions in this Style

## MathJax

MathJax is a JavaScript library that can render LaTeX math equations on a website. 
This template has MathJax  support enabled and can be used as shown below.

### Inline equation within the text

Sometimes we want to use a math expression within the text.  We can enclose
the LaTeX math expression within double-`$` characters like `$$ math-expression $$`.


This code:
~~~
Lorem ipsum dolor sit amet, consectetur  $$ k_{B}T/2 $$ adipisicing elit, sed 
do eiusmod tempor incididunt ut labore et dolore magna aliqua. 
~~~
{: .language-markdown}
... will be rendered like this:

Lorem ipsum dolor sit amet, consectetur  $$ k_{B}T/2 $$ adipisicing elit, sed 
do eiusmod tempor incididunt ut labore et dolore magna aliqua. 


### Math blocks

But we can also have our equations as blocks in between paragraphs. 

#### Centered Block Equation
This code:
~~~
$k_{n+1} = n^2 + k_n^2 - k_{n-1}$
{: .math-center}
~~~
{: .language-markdown}
... will align the equation to the center of the page:

$k_{n+1} = n^2 + k_n^2 - k_{n-1}$
{: .math-center}

#### Left Justified Block Equation
This code:
~~~
$k_{n+1} = n^2 + k_n^2 - k_{n-1}$
{: .math-left}
~~~
{: .language-markdown}
... will align the equation to the left of the page:

$k_{n+1} = n^2 + k_n^2 - k_{n-1}$
{: .math-left}


## Mermaid Examples

[Mermaid](https://mermaid-js.github.io/mermaid/#/) is a JavaScript library that can be used
to write various graphs and flow-charts within Markdown and render them on a website. 
This template has Mermaid support enabled and can be used as shown below.

### Flow Chart

```
<div class="mermaid">
graph LR
    A[Hard edge] -->|Link text| B(Round edge)
    B --> C{Decision}
    C -->|One| D[Result one]
    C -->|Two| E[Result two]
</div>
```

<div class="mermaid">
graph LR
    A[Hard edge] -->|Link text| B(Round edge)
    B --> C{Decision}
    C -->|One| D[Result one]
    C -->|Two| E[Result two]
</div>

### Git Graph

```
<div class="mermaid">
gitGraph:
options
{
    "nodeSpacing": 100,
    "nodeRadius": 10
}
end
commit
branch newbranch
checkout newbranch
commit
commit
checkout master
commit
commit
merge newbranch
</div>
```

<div class="mermaid">
gitGraph:
options
{
    "nodeSpacing": 100,
    "nodeRadius": 10
}
end
commit
branch newbranch
checkout newbranch
commit
commit
checkout master
commit
commit
merge newbranch
</div>
