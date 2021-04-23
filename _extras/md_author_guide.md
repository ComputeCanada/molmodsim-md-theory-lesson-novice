---
layout: page
title: "Authors Guide for Markdown Extensions"
permalink: /author-guide/
questions:
- "How can we add mathematical equations into our lessons?"
- "How can we add Flow diagrams as code within the lessons?"
objectives:
- "Demonstrate the use of MathJax equations inside a lesson."
- "Demonstrate the use of Mermaid charts inside a lesson."
keypoints:
- "Use `$$ LaTeX math $$` within a text."
- "Use `$ LaTeX math $` followed by `{: .math-center}` or `{: .math-left}` for math-blocks."
- 'Use a `<div class="mermaid"></div>` to add the code for a graph.'
---

* Table of Contents
{:toc}

## MathJax Examples

MathJax is a JavaScript library that can render LaTeX math equations on a website. 
This template has MathJax  support enabled and can be used as shown below.

### Inline equation within the text

Sometimes we want to use a math expression within the text.  We can enclose
the LaTeX math expression within double-`$` characters like `$$ math-expression $$`.

> ## Inline Math Expression
> This code:
> ~~~
> Lorem ipsum dolor sit amet, consectetur  $$ k_{B}T/2 $$ adipisicing elit, sed 
> do eiusmod tempor incididunt ut labore et dolore magna aliqua. 
> ~~~
> {: .language-markdown}
>
> ... will be rendered like this:
> 
> Lorem ipsum dolor sit amet, consectetur  $$ k_{B}T/2 $$ adipisicing elit, sed 
> do eiusmod tempor incididunt ut labore et dolore magna aliqua. 
>
{: .callout}


### Math blocks

But we can also have our equations as blocks in between paragraphs. 

> ## Centered Block Equation
> This code:
> ~~~
> $k_{n+1} = n^2 + k_n^2 - k_{n-1}$
> {: .math-center}
> ~~~
> {: .language-markdown}
>
> ... will align the equation to the center of the page:
> 
> $k_{n+1} = n^2 + k_n^2 - k_{n-1}$
> {: .math-center}
>
{: .callout}

> ## Left Justified Block Equation
> This code:
> ~~~
> $k_{n+1} = n^2 + k_n^2 - k_{n-1}$
> {: .math-left}
> ~~~
> {: .language-markdown}
>
> ... will align the equation to the left of the page:
> 
> $k_{n+1} = n^2 + k_n^2 - k_{n-1}$
> {: .math-left}
>
{: .callout}


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
