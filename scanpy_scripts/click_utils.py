"""
Provide helper functions for command line parsing with click
"""

import click


class CommaSeparatedText(click.ParamType):
    """
    Comma separated text
    """
    def __init__(self, dtype=click.STRING, simplify=False):
        self.dtype = dtype
        self.dtype_name = _get_type_name(dtype)
        self.name = '{}[,{}...]'.format(self.dtype_name, self.dtype_name)
        self.simplify = simplify

    def convert(self, value, param, ctx):
        try:
            converted = list(map(self.dtype, value.split(',')))
            if self.simplify and len(converted) == 1:
                converted = converted[0]
            return converted
        except ValueError:
            self.fail('{} is not a valid comma separated list of {}'.format(
                value, self.dtype_name), param, ctx)


def _get_type_name(obj):
    name = 'text'
    try:
        name = getattr(obj, 'name')
    except AttributeError:
        name = getattr(obj, '__name__')
    return name


def valid_limit(ctx, param, value):
    if value[0] > value[1]:
        param.type.fail(
            'lower limit must not exceed upper limit', param, ctx)
    return value


def valid_parameter_limits(ctx, param, value):
    for val in value:
        if val[1] > val[2]:
            param.type.fail(
                'lower limit must not exceed upper limit', param, ctx)
    return value


def mutually_exclusive_with(param_name):
    internal_name = param_name.strip('-').replace('-', '_').lower()
    def valid_mutually_exclusive(ctx, param, value):
        try:
            other_value = ctx.params[internal_name]
        except KeyError:
            return value
        if (value is None) == (other_value is None):
            param.type.fail(
                'mutually exclusive with "{}", one and only one must be '
                'specified.'.format(param_name),
                param,
                ctx,
            )
        return value
    return valid_mutually_exclusive


def required_by(param_name):
    internal_name = param_name.strip('-').replace('-', '_').lower()
    def required(ctx, param, value):
        try:
            other_value = ctx.params[internal_name]
        except KeyError:
            return value
        if other_value and not value:
            param.type.fail('required by "{}".'.format(param_name), param, ctx,)
        return value
    return required
