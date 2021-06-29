# hugo_api

## About

Scripts to handle the access to HUGO through their API that uses GraphQL

Reference: https://www.genenames.org/help/rest/

## Testing

This toolkit uses pytest and associated tests and data within the tests/ directory as part of its development.

These test should have close to full coverage across all the various objects, methods, and scripts.

    python3 -m pytest tests/data/test_hugo_utils.py

## Dependencies

### Python Modules
  requests
  requests_cache


# pharos_api

## About

Scripts to handle the access to Pharos through their API that uses GraphQL

Reference: https://pharos-api.ncats.io/graphql
GraphQL Reference: https://graphql.org/

## Testing

This toolkit uses pytest and associated tests and data within the tests/ directory as part of its development.

These test should have close to full coverage across all the various objects, methods, and scripts.

    python3 -m pytest tests/data/test_pharos_utils.py

## Dependencies

### Python Modules
  requests
  requests_cache
