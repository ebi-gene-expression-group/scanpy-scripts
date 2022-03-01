"""
Provide helper functions for command line parsing with click
"""

import click
import sys


class NaturalOrderGroup(click.Group):
    """Command group trying to list subcommands in the order they were added.

    With decorator, use::

        @click.group(cls=NaturalOrderGroup)
    """

    def list_commands(self, ctx):
        """List command names as they are in commands dict.

        If the dict is OrderedDict, it will preserve the order commands
        were added.
        """
        return self.commands.keys()


class CommaSeparatedText(click.ParamType):
    """
    Comma separated text
    """

    def __init__(self, dtype=click.STRING, simplify=False, length=None):
        self.dtype = dtype
        self.dtype_name = _get_type_name(dtype)
        self.simplify = simplify
        self.length = length
        if length and length <= 3:
            self.name = ",".join([f"{self.dtype_name}"] * length)
        else:
            self.name = "{}[,{}...]".format(self.dtype_name, self.dtype_name)

    def convert(self, value, param, ctx):
        """
        >>> @click.command()
        ... @click.option('--test-param')
        ... def test_cmd():
        ...     pass
        ...
        >>> ctx = click.Context(test_cmd)
        >>> param = test_cmd.params[0]
        >>> test_cst1 = CommaSeparatedText()
        >>> test_cst2 = CommaSeparatedText(click.INT, length=2)
        >>> test_cst3 = CommaSeparatedText(click.FLOAT, simplify=True)
        >>>
        >>> test_cst1.convert(None, param, ctx)
        >>> test_cst2.convert('7,2', param, ctx)
        [7, 2]
        >>> test_cst2.convert('7.2', param, ctx)
        Traceback (most recent call last):
        ...
        click.exceptions.BadParameter: 7.2 is not a valid integer
        >>> test_cst2.convert('7', param, ctx)
        Traceback (most recent call last):
        ...
        click.exceptions.BadParameter: 7 is not a valid comma separated list of length 2
        >>> test_cst3.convert('7.2', param, ctx)
        7.2
        """
        try:
            if value is None:
                converted = None
            else:
                converted = list(map(self.dtype, str(value).split(",")))
                if self.simplify and len(converted) == 1:
                    converted = converted[0]
        except ValueError:
            self.fail(
                "{} is not a valid comma separated list of {}".format(
                    value, self.dtype_name
                ),
                param,
                ctx,
            )
        if self.length:
            if len(converted) != self.length:
                self.fail(
                    "{} is not a valid comma separated list of length {}".format(
                        value, self.length
                    ),
                    param,
                    ctx,
                )
        return converted


class Dictionary(click.ParamType):
    """
    Text to be parsed as a python dict definition
    """

    def __init__(self, keys=None):
        self.name = "TEXT:VAL[,TEXT:VAL...]"
        self.keys = keys

    def convert(self, value, param, ctx):
        """
        >>> @click.command()
        ... @click.option('--my-param', type=Dictionary(keys=('abc', 'def', 'ghi', 'jkl', 'mno')))
        ... def test_cmd():
        ...     pass
        ...
        >>> ctx = click.Context(test_cmd)
        >>> param = test_cmd.params[0]
        >>> dict_param = param.type
        >>> dict_str1 = 'abc:0.1,def:TRUE,ghi:False,jkl:None,mno:some_string'
        >>> dict_str2 = 'abc:0.1,def:TRUE,ghi:False,jkl:None,mnp:some_string'
        >>> dict_str3 = ''
        >>> dict_param.convert(dict_str1, param, ctx)
        {'abc': 0.1, 'def': True, 'ghi': False, 'jkl': None, 'mno': 'some_string'}
        >>> dict_param.convert(dict_str2, param, ctx)
        Traceback (most recent call last):
        ...
        click.exceptions.BadParameter: mnp is not a valid key (('abc', 'def', 'ghi', 'jkl', 'mno'))
        >>> dict_param.convert(dict_str3, param, ctx)
        Traceback (most recent call last):
        ...
        click.exceptions.BadParameter:  is not a valid python dict definition
        """
        try:
            converted = dict()
            for token in value.split(","):
                if ":" not in token:
                    raise ValueError
                key, _, value = token.partition(":")
                if not key:
                    raise ValueError
                if isinstance(self.keys, (list, tuple)) and key not in self.keys:
                    self.fail(f"{key} is not a valid key ({self.keys})")
                if value == "None":
                    value = None
                elif value.lower() == "true":
                    value = True
                elif value.lower() == "false":
                    value = False
                else:
                    try:
                        value = float(value)
                    except ValueError:
                        pass
                converted[key] = value
            return converted
        except ValueError:
            self.fail(f"{value} is not a valid python dict definition", param, ctx)


def _get_type_name(obj):
    name = "text"
    try:
        name = getattr(obj, "name")
    except AttributeError:
        name = getattr(obj, "__name__")
    return name


def valid_limit(ctx, param, value):
    """
    Callback function that checks order of numeric inputs

    >>> @click.command()
    ... @click.option('--test-param', help='Sample help')
    ... def test_cmd():
    ...     pass
    ...
    >>> ctx = click.Context(test_cmd)
    >>> param = test_cmd.params[0]
    >>> valid_limit(ctx, param, value=[0.0125, 3])
    [0.0125, 3]
    >>> valid_limit(ctx, param, value=[0.0125, -0.0125])
    Traceback (most recent call last):
        ...
    click.exceptions.BadParameter: lower limit must not exceed upper limit
    >>> valid_limit(ctx, param, value=[0.0125, 0.0125])
    [0.0125, 0.0125]
    """
    if value[0] > value[1]:
        param.type.fail("lower limit must not exceed upper limit", param, ctx)
    return value


def valid_parameter_limits(ctx, param, value):
    """
    Callback function that checks order of multiple numeric inputs

    >>> @click.command()
    ... @click.option('--test-param', type=(click.STRING, click.FLOAT, click.FLOAT), multiple=True)
    ... def test_cmd():
    ...     pass
    ...
    >>> ctx = click.Context(test_cmd)
    >>> param = test_cmd.params[0]
    >>> valid_parameter_limits(ctx, param, [['a', 0.0, 2.0]])
    [['a', 0.0, 2.0]]
    >>> valid_parameter_limits(ctx, param, [['b', 0.0, 0.0]])
    [['b', 0.0, 0.0]]
    >>> valid_parameter_limits(ctx, param, [['c', 0.0, -1.0]])
    Traceback (most recent call last):
        ...
    click.exceptions.BadParameter: lower limit must not exceed upper limit
    >>> valid_parameter_limits(ctx, param, [['a', 0.0, 2.0], ['c', 0.0, -1.0]])
    Traceback (most recent call last):
        ...
    click.exceptions.BadParameter: lower limit must not exceed upper limit
    """
    for val in value:
        if val[1] > val[2]:
            param.type.fail("lower limit must not exceed upper limit", param, ctx)
    return value


def mutually_exclusive_with(param_name):
    internal_name = param_name.strip("-").replace("-", "_").lower()

    def valid_mutually_exclusive(ctx, param, value):
        try:
            other_value = ctx.params[internal_name]
        except KeyError:
            return value
        if (value is None) == (other_value is None):
            param.type.fail(
                'mutually exclusive with "{}", one and only one must be '
                "specified.".format(param_name),
                param,
                ctx,
            )
        return value

    return valid_mutually_exclusive


def required_by(param_name):
    internal_name = param_name.strip("-").replace("-", "_").lower()

    def required(ctx, param, value):
        try:
            other_value = ctx.params[internal_name]
        except KeyError:
            return value
        if other_value and not value:
            param.type.fail(
                'required by "{}".'.format(param_name),
                param,
                ctx,
            )
        return value

    return required


if __name__ == "__main__":
    import doctest

    sys.exit(doctest.testmod(verbose=True)[0])
