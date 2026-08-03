# Codex 结构化工具调用异常排障

这份文档给新会话使用，目标只有一个：当用户反馈“Codex 似乎会用 PowerShell，但命令被直接显示在聊天里，没有作为真实工具调用执行”时，先按这里排障，不要一上来就把问题归因到本机 PowerShell。

## 典型异常现象

下面这些都属于同一类问题：

- 助手消息里直接出现 `to=functions.exec`、`commentary to=functions.exec`、`to=functions.shell_command` 之类原始工具调用文本。
- 聊天中先出现一大段看起来像“工具调用 JSON / 指令体”的文字，然后下面才零散出现少量结果。
- 助手声称“我去执行 PowerShell”，但实际只是把命令文本发在对话框里。
- 某些旧会话或当前会话能正常调工具，但新建会话频繁把工具调用吐成纯文本。

这类问题的关键不是“命令有没有写对”，而是：**结构化工具调用在模型输出阶段就串格式了**。

## 先判断：是 shell 坏了，还是 tool-calling 串了

先区分两种情况。

### 情况 A：真实工具被调用了，但命令执行失败

特征：

- 能看到正常的工具执行结果块。
- 能拿到退出码、stdout、stderr。
- 只是命令本身报错，比如路径不存在、权限不足、PowerShell 语法错。

这种情况才是正常意义上的 shell / 命令问题。

### 情况 B：模型把工具调用文本直接发到聊天里

特征：

- 聊天正文里出现 `to=functions...`、`json`、`commentary to=...` 一类字面量。
- 看起来像“本该发给工具层的协议内容”直接显示给了用户。
- 有时同一条消息里既有原始工具调用文本，又混着部分执行结果。

这种情况优先判断为：

1. provider 与 Codex 的工具调用协议兼容不稳定；或
2. 当前模型 / 中转链路没有稳定返回结构化 tool call。

不要先把锅甩给 PowerShell。

## 本机已确认过的事实

以下结论来自这台机器的真实排查，可直接作为后续会话的初始判断依据：

1. 当前会话里 PowerShell 工具调用本身是可以正常工作的。
2. 异常新会话的本地 session 日志里，已经出现过字面量 `to=functions.exec` / `commentary to=functions.exec`。
3. 这说明至少有一部分问题不是命令执行层，而是**工具调用序列化/反序列化链路**。
4. 当前全局 Codex 配置使用的是第三方 provider：

```toml
model = "gpt-5.6-sol"
model_provider = "myproxy"
model_reasoning_effort = "low"

[model_providers.myproxy]
name = "MyProxy"
base_url = "https://sub2api.workhard.cyou/v1"
wire_api = "responses"
env_key = "MYPROXY_API_KEY"
```

这套组合的风险点不是单独某一项，而是三者叠加：

- 第三方代理 `myproxy`
- `wire_api = "responses"`
- 兼容性未充分验证的模型 slug：`gpt-5.6-sol`

## 最稳的判断顺序

用户一旦反馈这类问题，按下面顺序检查。

1. 先在当前会话里做一个极小的真实工具测试。
2. 再检查异常会话的本地 session jsonl 是否出现原始 `to=functions...` 文本。
3. 再看 `C:\Users\Administrator\.codex\config.toml` 的 provider / model / wire_api 组合。
4. 只有前面都排除了，才去怀疑 PowerShell 环境、路径、profile、系统权限。

## 最稳配置建议

下面分两档。

### A 档：最稳方案

如果用户最看重的是“Codex 稳定调用 PowerShell / 文件工具 / 浏览器工具”，优先建议：

- 不要走第三方 proxy 作为 Codex 的主 provider。
- 暂时移除自定义 `myproxy` provider 配置。
- 使用 Codex Desktop 官方登录/默认 provider。

原因很直接：

**要降低工具调用异常，最有效的办法就是减少兼容层。**

也就是：

- 官方 provider > 第三方 proxy
- 已验证稳定的主流模型 > 兼容性不透明的自定义/新 slug

### B 档：必须保留代理时的保守方案

如果用户必须继续使用代理：

- 保留 `model_provider = "myproxy"`
- 但不要默认继续使用 `gpt-5.6-sol`
- 改成代理文档里明确支持、且已经实测能稳定 tool-calling 的标准模型 ID
- 只有确认代理正确实现了 Codex 所需结构化工具调用，才继续保留 `wire_api = "responses"`

这里要明确告诉用户：

**如果代理本身对 Responses API 的结构化工具调用支持不完整，本地 `config.toml` 再怎么改，也只能缓解，不能根治。**

## 不要优先折腾的项

除非有额外证据，否则不要把主要精力先放在这些地方：

- PowerShell 本体重装
- 改 `sandbox = "elevated"`
- 关闭插件
- 怀疑项目目录权限
- 怀疑工作区内容导致工具完全不可用

这些都不是这次问题里最强的证据方向。

## 推荐给用户的验证步骤

改完配置后，让用户这样验证：

1. 完全退出 Codex Desktop。
2. 重新打开后，新建一个全新线程。
3. 发送：`请直接执行 Get-Location 并返回结果`
4. 再发送：`请执行 Write-Output 'ok'; Get-Date -Format o`
5. 观察聊天正文里是否还出现裸露的 `to=functions...` 文本。

如果仍然出现，大概率还是 provider/tool-call 链路不稳。

## 给新会话的结论模板

当证据已经足够时，可以直接这样总结：

> 这更像 provider/tool-calling 协议兼容问题，不像本机 PowerShell 损坏。当前证据显示，异常会话把原本应该发给工具层的结构化调用直接渲染成了聊天文本。最稳的办法是先切回官方 provider；如果必须保留代理，再换成代理已验证稳定的标准模型 ID，并重新做新线程测试。

## 这份文档的边界

这份文档只解决“为什么工具调用被显示成聊天文本”这一类问题，不负责判断具体 PowerShell 命令语义是否正确，也不替代项目内的业务排障文档。
