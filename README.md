# BootsCrispr



## Introduction
BootsCrispr is a web-based tool used for sgRNA assisted design. In addition to assisting in sgRNA assisted design, it can
also search for potential off target locations within sequences.

The current version is compatible with three network models: CNN, RNN, and Transformer.

 
 
## Requirement
* django=3.2.15
* python == 3.7
* tensorflow == 1.13.1
* sonnet == 1.9
* keras=2.2.4
* Cas-OFFFinder = 3.0


Note:
1. BootCrispr currently provides visual pages for users to use ;
2. In addition to conventional species, users can provide us with species information that needs to be initialized for user data analysis and design.

## Usage
1. home

![home](static/images/img_6.png)

2. sgRNA Design

![design1](static/images/Figure5.png)
![design2](static/images/Figure6-A.png)
![design3](static/images/Figure6-B.png)
![design4](static/images/Figure7.png)

3. Search Off Targets

![design5](static/images/Figure8.png)
![design6](static/images/Figure9-A.png)
![design7](static/images/Figure9-B.png)


4. About
![about](static/images/about.jpg)


#### Prediction
There are three models: RNN, Transformer, and off target models, including CNN network models.

## Citation

## Contacts
bioinfo2025@163.com
