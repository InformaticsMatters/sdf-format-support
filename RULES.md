# Rules

>   In this text the terms **must**, **must not** and **should** should be
    interpreted according to best practices, as described in [RFC-2119].

To properly comply with our requirements and assist in future
automation techniques you **must** satisfy the following design rules.

>   In order to provide a formatter that is 'approved' by Informatics Matters
    we may introduce automated background testing and inspection of your
    repository files, and the images you publish by inspecting this repository,
    your test data and Docker Hub. We may also introduce our own automated
    build of your images.

>   Deviating from the rules laid out in this file may result in your formatter
    being removed from the Data Tier Manager format library.

In summary, for a format-support container image to run successfully
within the Data Tier Management service...

-   It **must** have a `:stable` image tag
-   It **must not** exit with a non-zero exit code
-   It **should** write something to `DT_DATASET_OUTPUT_PATH`

## Repository design rules

Repositories...

1.  **Must** not remove or edit the files `RULES.md`, `TESTING.md`,
    `EVENTS.md`, `.yamllint` or `VERSION.txt`, these files are excluded from
    the MIT license and must remain, unaltered if the formatter is to pass
    potential future automated testing
2.  **Must** not remove the GitHub Action `lint` Job in any GitHub
    Action workflow, i.e. the `lint` action is an acceptance requirement
3.  **Must** be [created] using this template repository 
4.  **Must** root any and all source code (implementation) in
    the `source` directory
5.  **Should** place format support documentation in the project `README.md`
6.  **Must** be named `<type>-format-support`, where _type_ is a
    symbolic reference of the type of dataset you're supporting (i.e `sdf`)
7.  **Must** produce public images on Docker Hub that can be accessed using
    the reference `<owner>/<type>-format-support:<tag>`
8.  **Must** be buildable using a `Dockerfile` in the root of the project
9.  **Must** support the production of an image using the tag `stable`
10. **Should** support the production of image tags for `latest`
    and individual [Semantic Versioning 2] tags like `1.0.0`
11. **Must** support being built using the command `docker build .`
    from the project root, without additional parameters
12. **Must** contain at least one test dataset that can be processed
    successfully
13. **Must** store test data that can be processed successfully in
    a sub-directory of `test/success`, i.e. `test/success/1/good.data`
14. **Should** contain at least one test dataset file that cannot be processed
    successfully (i.e. an error-handling test)
15. **Must** store test data that cannot be processed successfully in
    a sub-directory of `test/failure`, i.e. `test/failure/1/bad.data`
16. **Must** provide `build`, `build latest`, `build tag` and `build stable`
    badges on a single line, in format provided in the top of the current
    `README.md`
17. **Must** provide a `tag` badge, on its own line, in format provided
    in the top of the current `README.md`
18. **Should** have access to docker-compose v1.27.0 or better
19. **Should** have access to Python 3.8 or better - especially if the
    developer wants to execute the GitLab Action lint tests
 
>   Read TESTING.md for a discussion of the test directory structure
    and example test execution commands that you are expected to be
    able to support

>   Read EVENTS.md for a discussion of the event-reporting mechanism
    that can be used to pass information back to the user.
 
## Image behaviour rules

Images...

1.  **Must** be published as a public image to Docker Hub
2.  **Should** use `docker-entrypoint.sh` as their container entrypoint
3.  **Must** place the implementation into the directory `/home/format-support`
    as the Pod will be started with this as the working directory.
4.  **Must** use `CMD` and not `ENTRYPOINT` to launch the entrypoint script  
5.  **Must** expect to be executed using an arbitrary user and group ID.
    Your container cannot expect to run as a privileged user
6.  **Must** expect the environment variable `DT_DATASET_NAME` to be set
    to a string representing the name given of the dataset provided
7.  **Must** expect the environment variable `DT_DATASET_FILE` to be set
    to a string representing the full path to the dataset file that is to be
    processed
8.  **Must** process the input dataset into files in the directory
    identified by `DT_DATASET_OUTPUT_PATH`
9.  **Should** expect to be limited to no more than 1 CPU core ane no
    more than 1GiB of memory. Importantly, exceeding the memory limit will
    result in the container being terminated
10. **Must** use a non-zero exit code to indicate an unrecoverable failure

---

[created]: https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/creating-a-repository-from-a-template
[rfc-2119]: https://tools.ietf.org/html/rfc2119
[semantic versioning 2]: https://semver.org
