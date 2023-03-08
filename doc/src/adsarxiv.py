"""
ADSArxiv PybTeX style
"""
# Copyright Emmanuel Bertin CFHT/CNRS/SorbonneU
# Licensed under GPL v3

from packaging import version as vers
from pybtex import __version__ as pybtex_version
from pybtex.style.formatting.unsrt import Style as UnsrtStyle, date, pages, toplevel
from pybtex.style.template import * # ... and anything else needed
from pybtex.plugin import register_plugin

pybtex_new_version = vers.parse(pybtex_version) >= vers.parse("0.22")

if pybtex_new_version:
    def _format_data(node, data):
        try:
            f = node.format_data
        except AttributeError:
            return node
        else:
            return f(data)


    def _format_list(list_, data):
        return (_format_data(part, data) for part in list_)

    @node
    def href2(children, data):
        parts = _format_list(children, data)
        if "http" in list(parts)[0]:
            parts = _format_list(children, data)
            return richtext.HRef(*parts)
        else:
            parts = _format_list(children, data)
            return richtext.Tag('strong', *parts)

class ADSArxivStyle(UnsrtStyle):

    if pybtex_new_version:
        def format_names(self, role, as_sentence=True):
            formatted_names = names(role, sep=', ', sep2 = ' and ', last_sep=', and ')
            if as_sentence:
                return sentence[ formatted_names ]
            else:
               return formatted_names

        def get_article_template(self, e):
            volume_and_pages = first_of[
            # volume and pages, with optional issue number
                optional[
                    join[
                        field('volume'),
                        optional[ '(', field('number'),')' ],
                        ':', pages
                    ],
                ],
                # pages only
                words[ 'pages', pages ],
            ]
            myurl = first_of[
                    optional_field('adsurl'),
                    optional[ join[ 'http://arxiv.org/abs/', field('eprint') ]],
                    optional_field('url'),
                    optional[ join ['http://dx.doi.org/', field('doi') ]]
            ]
            template = toplevel[
                sentence[ self.format_names('author', as_sentence=False), field('year') ],
                href2[ myurl, self.format_title(e, 'title') ],
                sentence(capfirst=False) [
                    tag('emph')[ field('journal') ],
                    optional[ volume_and_pages ]],
                sentence(capfirst=False) [ optional_field('note') ],
            ]
            return template

        def get_book_template(self, e):
            myurl = first_of[
                optional_field('adsurl'),
                optional[ join ['http://arxiv.org/abs/', field('eprint') ]],
                optional_field('url'),
                optional[ join ['http://dx.doi.org/', field('doi') ]]
            ]
            template = toplevel[
                self.format_author_or_editor(e),
                href2[ myurl, self.format_btitle(e, 'title') ],
                self.format_volume_and_series(e),
                sentence [
                    field('publisher'),
                    optional_field('address'),
                    self.format_edition(e),
                    date
                ],
                optional[ sentence [ self.format_isbn(e) ] ],
                sentence [ optional_field('note') ],
            ]
            return template

        def get_inproceedings_template(self, e):
            myurl = first_of[
                optional_field('adsurl'),
                optional[ join ['http://arxiv.org/abs/', field('eprint') ]],
                optional_field('url'),
                optional[ join ['http://dx.doi.org/', field('doi') ]]
            ]
            template = toplevel[
                sentence[ self.format_names('author', as_sentence=False), field('year') ],
                href2[ myurl, self.format_title(e, 'title') ],
                words[
                    'In',
                    sentence(capfirst=False)[
                        optional[ self.format_editor(e, as_sentence=False) ],
                        self.format_btitle(e, 'booktitle', as_sentence=False),
                        self.format_volume_and_series(e, as_sentence=False),
                        optional[ pages ],
                    ],
                    self.format_address_organization_publisher_date(e),
                ],
                sentence(capfirst=False)[ optional_field('note') ],
            ]
            return template

        def get_misc_template(self, e):
            myurl = first_of[
                    optional_field('adsurl'),
                    optional[ join ['http://arxiv.org/abs/', field('eprint') ]],
                    optional_field('url'),
                    optional[ join ['http://dx.doi.org/', field('doi') ]]
            ]
            template = toplevel[
                optional[ sentence[ self.format_names('author', as_sentence=False), optional [ field('year') ]]],
                optional[ href2[ myurl, self.format_title(e, 'title') ]],
                sentence[ optional[ field('howpublished') ]],
                sentence[ optional_field('note') ],
            ]
            return template
    else:

        def format_article(self, e):
            volume_and_pages = first_of [
            # volume and pages, with optional issue number
                optional [
                    join [
                        field('volume'),
                        optional['(', field('number'),')'],
                        ':', pages
                    ],
                ],
                # pages only
                words ['pages', pages],
            ]
            myurl = first_of [
                    optional_field('adsurl'),
                    optional [join ['http://arxiv.org/abs/'], field('eprint')],
                    optional_field('url'),
                    optional [join ['http://dx.doi.org/', field('doi')]]
            ]
            template = toplevel [
                self.format_names('author'),
                href [myurl, self.format_title(e, 'title')],
                sentence(capfirst=False) [
                    tag('emph') [field('journal')],
                    optional[ volume_and_pages ],
                    field('year')],
                sentence(capfirst=False) [ optional_field('note') ],
            ]
            return template.format_data(e)

        def format_book(self, e):
            myurl = first_of [
                optional_field('adsurl'),
                optional [join ['http://arxiv.org/abs/'], field('eprint')],
                optional_field('url'),
                optional [join ['http://dx.doi.org/', field('doi')]]
            ]
            template = toplevel [
                self.format_author_or_editor(e),
                href [myurl, self.format_btitle(e, 'title')] \
                    if len(myurl.format_data(e)) > 0 \
                    else tag('strong') [self.format_btitle(e, 'title')],
                self.format_volume_and_series(e),
                sentence [
                    field('publisher'),
                    optional_field('address'),
                    self.format_edition(e),
                    date
                ],
                optional[ sentence [ self.format_isbn(e) ] ],
                optional_field('note'),
            ]
            return template.format_data(e)

        def format_inproceedings(self, e):
            myurl = first_of [
                optional_field('adsurl'),
                optional [join ['http://arxiv.org/abs/', field('eprint')]],
                optional_field('url'),
                optional [join ['http://dx.doi.org/', field('doi')]]
            ]
            template = toplevel [
                sentence [self.format_names('author')],
                href [myurl, self.format_title(e, 'title')] \
                    if len(myurl.format_data(e)) > 0 \
                    else tag('strong') [self.format_title(e, 'title')],
                words [
                    'In',
                    sentence(capfirst=False) [
                        optional[ self.format_editor(e, as_sentence=False) ],
                        self.format_btitle(e, 'booktitle', as_sentence=False),
                        self.format_volume_and_series(e, as_sentence=False),
                        optional[ pages ],
                    ],
                    self.format_address_organization_publisher_date(e),
                ],
                sentence(capfirst=False) [ optional_field('note') ],
            ]
            return template.format_data(e)


register_plugin('pybtex.style.formatting', 'adsarxiv', ADSArxivStyle)

