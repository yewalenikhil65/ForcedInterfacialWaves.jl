/* Pair the existing Julia and MATLAB code blocks on theory pages into a single
 * tabbed widget: readers click a "Julia" or "MATLAB" button to reveal the
 * corresponding code (only one is shown at a time).
 *
 * Executable Julia examples may generate a <pre> output block or an image
 * before the MATLAB block; those nodes stay with the Julia panel.
 */
(function () {
    function isLanguageBlock(node, language) {
        return node && node.matches("pre") &&
            node.querySelector("code.language-" + language);
    }

    function isGeneratedJuliaOutput(node) {
        if (!node) {
            return false;
        }
        // Executed @example output, or images/svg produced by a plot.
        if (node.matches("pre.documenter-example-output") ||
            node.matches("img") ||
            node.matches("svg") ||
            node.classList.contains("documenter-example-output")) {
            return true;
        }
        // A plain ```text``` (or other non-matlab) fenced block placed between the
        // Julia and MATLAB code to show printed output. Documenter renders it as
        // <pre><code class="language-text hljs">. Keep it with the Julia panel.
        if (node.matches("pre")) {
            var code = node.querySelector("code");
            if (!code) {
                return true;
            }
            if (!code.classList.contains("language-julia") &&
                !code.classList.contains("language-matlab")) {
                return true;
            }
        }
        return false;
    }

    function isBoundary(node) {
        return node && node.matches("h1, h2, h3, h4, h5, h6, hr, p, ul, ol, table, figure");
    }

    function makeTabButton(label, panel, group) {
        var button = document.createElement("button");
        button.type = "button";
        button.className = "code-tab-button";
        button.textContent = label;
        button.addEventListener("click", function () {
            group.buttons.forEach(function (b) { b.classList.remove("active"); });
            group.panels.forEach(function (p) {
                p.classList.remove("active");
                p.style.display = "none";
            });
            button.classList.add("active");
            panel.classList.add("active");
            panel.style.display = "block";
        });
        return button;
    }

    function pairExamples() {
        document.querySelectorAll("pre").forEach(function (juliaBlock) {
            if (!isLanguageBlock(juliaBlock, "julia") ||
                juliaBlock.closest(".code-comparison")) {
                return;
            }

            var matlabBlock = null;
            var generatedOutput = [];
            var node = juliaBlock.nextElementSibling;

            while (node) {
                if (isLanguageBlock(node, "matlab")) {
                    matlabBlock = node;
                    break;
                }
                if (isBoundary(node) || !isGeneratedJuliaOutput(node)) {
                    return;
                }
                generatedOutput.push(node);
                node = node.nextElementSibling;
            }

            if (!matlabBlock) {
                return;
            }

            // Collect any generated-output blocks (e.g. a ```text``` printout)
            // that immediately follow the MATLAB block, so the MATLAB tab shows
            // its own output rather than leaving it always-visible below the tabs.
            var matlabOutput = [];
            var after = matlabBlock.nextElementSibling;
            while (after && !isBoundary(after) && isGeneratedJuliaOutput(after)) {
                matlabOutput.push(after);
                after = after.nextElementSibling;
            }

            // Build the tab container.
            var wrapper = document.createElement("div");
            wrapper.className = "code-comparison code-tabs";

            var tabBar = document.createElement("div");
            tabBar.className = "code-tab-bar";
            tabBar.setAttribute("role", "tablist");

            var juliaPanel = document.createElement("div");
            juliaPanel.className = "code-tab-panel code-tab-panel-julia active";
            juliaPanel.style.display = "block";

            var matlabPanel = document.createElement("div");
            matlabPanel.className = "code-tab-panel code-tab-panel-matlab";
            matlabPanel.style.display = "none";

            var group = { buttons: [], panels: [juliaPanel, matlabPanel] };

            var juliaButton = makeTabButton("Julia", juliaPanel, group);
            juliaButton.classList.add("code-tab-julia", "active");
            var matlabButton = makeTabButton("MATLAB", matlabPanel, group);
            matlabButton.classList.add("code-tab-matlab");
            group.buttons = [juliaButton, matlabButton];

            tabBar.appendChild(juliaButton);
            tabBar.appendChild(matlabButton);

            // Move the blocks into their panels.
            juliaBlock.parentNode.insertBefore(wrapper, juliaBlock);
            wrapper.appendChild(tabBar);
            wrapper.appendChild(juliaPanel);
            wrapper.appendChild(matlabPanel);

            juliaPanel.appendChild(juliaBlock);
            generatedOutput.forEach(function (output) {
                juliaPanel.appendChild(output);
            });
            matlabPanel.appendChild(matlabBlock);
            matlabOutput.forEach(function (output) {
                matlabPanel.appendChild(output);
            });
        });
    }

    // Documenter bundles highlight.js (11.8.0) with the Julia and julia-repl
    // grammars, but not MATLAB, so MATLAB code fences render uncoloured. Load the
    // matching MATLAB grammar from the same CDN/version, then (re)highlight every
    // MATLAB block so it gets the same token colouring as the Julia blocks.
    var HLJS_VERSION = "11.8.0";
    var MATLAB_GRAMMAR_URL =
        "https://cdnjs.cloudflare.com/ajax/libs/highlight.js/" +
        HLJS_VERSION + "/languages/matlab.min.js";

    function highlightMatlabBlocks() {
        if (typeof hljs === "undefined") {
            return;
        }
        document.querySelectorAll("pre code.language-matlab").forEach(function (block) {
            if (block.dataset.highlighted === "yes") {
                return;
            }
            hljs.highlightElement(block);
        });
    }

    function enableMatlabHighlighting() {
        // hljs must be present; poll briefly because Documenter loads it via
        // require.js and it may not be ready the instant this script runs.
        var attempts = 0;
        var timer = setInterval(function () {
            attempts += 1;
            if (typeof hljs !== "undefined") {
                clearInterval(timer);
                if (!hljs.getLanguage("matlab")) {
                    var script = document.createElement("script");
                    script.src = MATLAB_GRAMMAR_URL;
                    script.onload = highlightMatlabBlocks;
                    document.head.appendChild(script);
                } else {
                    highlightMatlabBlocks();
                }
            } else if (attempts > 40) {
                clearInterval(timer); // give up after ~10 s
            }
        }, 250);
    }

    function run() {
        pairExamples();
        enableMatlabHighlighting();
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", run);
    } else {
        run();
    }
})();
