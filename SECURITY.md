# Security policy

## Supported version

Only the current `main` branch is maintained. This is an educational simulation
and benchmark study, not a production astrophysics service.

## Reporting a vulnerability

Please use GitHub's private vulnerability-reporting form under the repository's
Security tab. Do not open a public issue for a report that contains an exploit,
credential, or other sensitive detail. Ordinary non-sensitive defects can be
reported through GitHub Issues.

The executable reads command-line arguments and optionally writes snapshot CSV
files. Do not run untrusted builds or supply output paths that you do not intend
the process to create.
