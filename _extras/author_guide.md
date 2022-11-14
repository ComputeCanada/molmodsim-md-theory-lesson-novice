---
layout: page
title: "Authors Guide"
permalink: /author-guide/
questions:
- "How can I create lessons using this style?"
- "How do I update the style of my lesson?"
- "How can we add mathematical equations into our lessons?"
- "How can we add Flow diagrams as code within the lessons?"
- "How can I hide blocks of text and instead show a summary during a workshop?"
objectives:
- "Give a quick overview on how this style is used."
- "Demonstrate the use of MathJax equations inside a lesson."
- "Demonstrate the use of Mermaid charts inside a lesson."
- "Demonstrate _Instructor View_."
keypoints:
- "Use `$$ LaTeX math $$` within a text."
- "Use `$ LaTeX math $` followed by `{: .math-center}` or `{: .math-left}` for math-blocks."
- 'Use a `<div class="mermaid"></div>` to add the code for a graph.'
- "Mark paragraphs or blocks with `{: .self_study_text :}` or  `{: .instructor_notes :}`."
---

* Table of Contents
{:toc}

# Preface

[This style]({{ site.cc_style_repo }}) is based on based on The Carpentries style and has some 
modifications for use by the [Digital Research Alliance of Canada]({{ site.alliance_site }}) 
(the Alliance) and [ACENET]({{ site.an_site }}). 
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
* `carpentry`:  This changes the branding of the lesson, e.g. `carpentry: "alliance"`
   for the Alliance.
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

### Flowchart

[Documentation for Mermaid Flowchart](https://mermaid-js.github.io/mermaid/#/flowchart)

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

[Documentation for Mermaid gitGraph](https://mermaid-js.github.io/mermaid/#/gitgraph)

```
<div class="mermaid">
gitGraph
   commit id: "C1"
   commit id: "C2"
   branch develop
   checkout develop
   commit id: "C3"
   commit id: "C4"
   checkout main
   merge develop
   commit id: "C5"
   commit id: "C6"
</div>
```

<div class="mermaid">
gitGraph
   commit id: "C1"
   commit id: "C2"
   branch develop
   checkout develop
   commit id: "C3"
   commit id: "C4"
   checkout main
   merge develop
   commit id: "C5"
   commit id: "C6"
</div>

## Instructor View

Some lessons are rather complex and require a lot of information, which can result in long blocks
of text that are needed for self-study as well as new instructors that need to familiarize
themselves with the material.  
In a workshop or presentation setting however it's difficult to pick out the talking points and
having key information in bullet-points would be much better.

A solution to this is to introduce an "Instructor View" that will hide paragraphs that have been 
marked as `{: .self_study_text :}` and instead show paragraphs marked as `{: .instructor_notes :}`,
which are otherwise hidden.  

If any of those markers a present on a particular page (and only on those), the following toggle
switch will appear in the Navigation bar at the top to turn the _instructor view_ on or off:

<div class="toggle_instructor_view" title="Toggle between self-study text and instructor notes" data-toggle="tooltip">
  <span class="glyphicon glyphicon-education" aria-hidden="true"></span>&nbsp;<label class="switch">
    <input type="checkbox">
    <span class="slider round"></span>
  </label>
</div>

The state of the toggle is persistent for some time (by default 7 days) by setting a cookie, 
however this duration can be customized by a variable in `_config.yml`:

```yaml
instructor_view_cookie_lifetime: 3
```


**Example**:

```
This is a long text that describes a complex topic in detail. 
This contains information that the instructor explains in detail.   
_Lorem ipsum dolor sit amet,[...] sunt in culpa qui officia deserunt mollit anim id est laborum._
{: .self_study_text :}

* talking points for instructor
* easier to read during the workshop
{: .instructor_notes :}
```

This is a long text that describes a complex topic in detail. 
This contains information that the instructor explains verbally in detail.   
_Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut_
_labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris_
_nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate_
_velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non_
_proident, sunt in culpa qui officia deserunt mollit anim id est laborum._
_Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut_
_labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris_
_nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate_
_velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non_
_proident, sunt in culpa qui officia deserunt mollit anim id est laborum._
{: .self_study_text :}

* talking points for instructor
* easier to read during the workshop
{: .instructor_notes :}

