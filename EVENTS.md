# Format support event reporting
Format containers can report status messages back to the data tier as _events_,
which will make them available via the REST API. Status messages that you want
to be recorded as events should be written to the container's stdout.

>   As a general guide, developers should otherwise limit what's written to
    stdout from a format image - as each line is processed by the DataTier.

## Event line format
Messages that you want to be interpreted as _events_ are lines in stdout
with the following format: -

    <timstamp> # <event type> -EVENT- <message>

Where: -

-   `timestamp` is any string representing a date and time that can be
    processed by v2.8.x of the Python [python-dateutil] module
-   `event type` is one of a standard set of recognised Python logging [levels]
-   `message` is any short message, which will be truncated to 79 characters

Lines found in stdout that do not match this format will be ignored.

## Interpretation of event types
At the moment the event types (logging levels) are all interpreted
in the same way. The event type is preserved and returned to the user
along with its message. Importantly, presence of error-like types (`CRITICAL`
and `ERROR` for example) are not interpreted as a failure of the format-support
container. Regardless of the presence of a `CRITICAL` event type, the container
is only considered to have failed if it uses a non-zero exit code.

## Examples
The following lines all represent valid event lines: -

-   `2021-02-25T12:23:00.123456Z # INFO -EVENT- 400 of 25,000 molecules`
-   `Sat Oct 11 17:13:46 UTC 2003 # WARNING -EVENT- No molecules`

---

[python-dateutil]: https://pypi.org/project/python-dateutil/
[levels]: https://docs.python.org/3/library/logging.html#logging-levels
