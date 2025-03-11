# MitoQuest Changelog

## [TBD] - TBD

- RAII 原则
- 重写并完善了 iobgzf.h 中的 BGZFile class，实现对 htslib/bgzf.h 的封装
- 使用 BGZFile class 不再直接调用 BGZF，这样更安全、更方便，也符合 RAII 原则

## [1.3.0] - 2025-03-10

- 创建 Github Actions workflows 实现跨平台编译
- 发布 Release 时自动附带编译好的二进制可执行文件
