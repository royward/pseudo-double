# Security Policy

Please don't use this code anywhere mission critical - it has not been tested anywhere near sufficiently for that.

This code does no heap allocation, uses no pointers, uses no libraries and makes no system calls. The only potential security implications are if something was caused downstream by the results not being as accurate as expected.

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 1.0.2   | :white_check_mark: |
| 1.0.1   | :x:                |