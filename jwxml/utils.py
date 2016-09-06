__all__ = ['iterchildren', 'etree']

try:
    from lxml import etree

    def iterchildren(element, tag=None):
        return element.iterchildren(tag)
except ImportError:
    import xml.etree.cElementTree as etree

    # The ElementTree implementation in xml.etree does not support
    # Element.iterchildren, so provide this wrapper instead
    # This wrapper does not currently provide full support for all the arguments as
    # lxml's iterchildren
    def iterchildren(element, tag=None):
        if tag is None:
            return iter(element)

        def _iterchildren():
            for child in element:
                if child.tag == tag:
                    yield child

        return _iterchildren()
