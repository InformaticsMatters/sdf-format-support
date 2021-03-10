# format-support-template

[![build](https://github.com/InformaticsMatters/format-support-template/actions/workflows/build.yaml/badge.svg)](https://github.com/InformaticsMatters/format-support-template/actions/workflows/build.yaml)
[![publish latest](https://github.com/InformaticsMatters/format-support-template/actions/workflows/publish-latest.yaml/badge.svg)](https://github.com/InformaticsMatters/format-support-template/actions/workflows/publish-latest.yaml)
[![publish tag](https://github.com/InformaticsMatters/format-support-template/actions/workflows/publish-tag.yaml/badge.svg)](https://github.com/InformaticsMatters/format-support-template/actions/workflows/publish-tag.yaml)
[![publish stable](https://github.com/InformaticsMatters/format-support-template/actions/workflows/publish-stable.yaml/badge.svg)](https://github.com/InformaticsMatters/format-support-template/actions/workflows/publish-stable.yaml)

![GitHub tag (latest SemVer)](https://img.shields.io/github/v/tag/informaticsmatters/format-support-template)

---

>   Replace the following template text with your own documentation

A template repository for **Informatics Matters DataTier dataset
format-support container image implementations**.

Use this repository to start development of a container image that can be
used by the DataTier Manager to support the loading of a dataset _type_.

You are required to preserve and adhere to the design rules
that can be found in this repository's `RULES.md`. 

Please take time to familiarise yourself with following sections before
discarding this text. Failing to comply with what's discussed here
may result in the image not passing any future automated acceptance testing
that may be deployed.

-   Building
-   Testing
-   Required image tags
-   Built-in GitHub Actions
-   Updating the badge links

## Building
You should be able to build your format-support image using `docker`: -

    $ docker build . -t ${PWD##*/}:latest

Or `docker-compose`: -

    $ IMAGE_NAME=${PWD##*/} docker-compose build

>   The expectation is that you're using a recent Docker engine
    (like `20.10.2`) and docker-compose utility (like `1.27.4`).

## Testing
Familiarise yourself with `TESTING.md`, which describes
a strategy you can replicate while developing and testing your image.

## Required image tags
The DataTier Manager will only execute formatter images tagged `:stable`
so you **MUST** ultimately produce an image with this tag. You are also
encouraged to produce a `:latest` tag and any specific tags that satisfy your
own needs (with formats like `1.0.0-rc.1`, `1.0.0` and `2021.1`).

>   This repository's built-in GitHub Actions (see the next section)
    will do all this for you.

## Built-in GitHub Actions
The template contains a number of [GitHub Actions] that will automatically
build the container image and also publish `:latest` and any tags you make to
Docker Hub.

This relies on your docker registry mirroring your repository. If your docker
repository name does not mirror your GitHub repository name then you will need
to adjust these actions. For example, if you create a GitHub repository
from this one and call it 'XYZ/my-support-template' then you must be able to
push docker images to 'xyz/my-support-template:latest'. If not, you will need
to edit the workflow files to satisfy your needs.

>   Your images must be published to Docker Hub.

The following built-in actions are: -

1.  **For every commit to main** an Action builds the docker image and
    pushes it using the image tag `:latest`.
    This is accomplished by the `publish-latest.yaml` Action.
2.  **For every repository tag** an Action builds the docker image and
    pushes it with the image tag `:<tag>`
    This is accomplished by the `publish-tag.yaml` Action.
    -   If the tag looks like a _formal_ release, i.e. is a 2 or 3-digit number
        like `2021.1` or `1.0.0` an Action builds the docker image and
        pushes it using the image tag `:<tag>` and the tag `:stable`
        This is accomplished by the `publish-stable.yaml` Action.
3.  **For every commit on a branch**, or a pull request to main, an Action
    runs that just builds the docker image - but does not push it.
    This is accomplished by the `build.yaml` Action.

In order for the above Actions to succeed you will need to define the following
GitHub Repository (or Organisation) **secrets**: -

-   `DOCKERHUB_USERNAME` A Docker Hub user
-   `DOCKERHUB_TOKEN` A Docker Hub user password or, ideally, access token

>   Repositories created in the InformaticsMatters organisation
    are already presented with these secrets as they are already
    defined at the Organisation level.

## Updating the badge links
Don't forget to replace the `InformaticsMatters/format-support-template`
values in the above badge links with the name of your own repository, otherwise
your badges will reflect the template repository's state, not yours.

---

[github actions]: https://github.com/features/actions
