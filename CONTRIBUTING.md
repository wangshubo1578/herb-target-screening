# 贡献指南 / Contributing

欢迎对本项目提出改进！中文说明如下（English below）。

## 如何贡献

1. **Fork** 本仓库到你自己的账号
2. 创建你的分支：`git checkout -b feature/your-feature`
3. 提交修改：`git commit -m "Add: 你的改动说明"`
4. 推送分支：`git push origin feature/your-feature`
5. 发起 **Pull Request**

## 贡献方向建议

- 新增中药数据源（如 SymMap、BATMAN-TCM 等的接入）
- 优化化合物清洗规则（`assets/filter_blacklist.json`、`compound_aliases.json`）
- 改进对接协议（盒子定义、打分函数参数）
- 报告可视化增强
- 修复 bug、补充文档

## 代码规范

- Python 脚本保持单文件自包含风格，遵循 PEP 8
- 每个步骤脚本保持「输入 JSON → 输出 JSON」的接口约定（见 SKILL.md 流程图）
- 新增依赖需同步更新 `requirements.txt` 与 README

---

## English

Contributions are welcome! Please fork the repository, create a feature
branch, and open a Pull Request. Keep step scripts self-contained with the
JSON-in / JSON-out interface convention described in SKILL.md. Update
`requirements.txt` and the README when adding dependencies.
