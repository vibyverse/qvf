Prepare the documentation site for publication (WP-64). `mkdocs build
--strict` with Doxygen is clean: the Doxyfile no longer depends on Graphviz
(`HAVE_DOT = NO`), a `tensor_view.hpp` comment no longer names a zero-argument
`Writer::add()` that does not exist, and strict mode now also checks heading
anchors. A dispatch-only GitHub Pages workflow (`docs-pages.yml`, pinned
actions) builds and deploys the site. CONTRIBUTING.md lists the one-time
Pages settings. Publication stays a maintainer action and has not happened.
